# Phase 2c: Direct Modification Calling by Sample Comparison
# Enabled by '^mod-calls <backend>' in CONFIG (uncalled4 and/or remora).
#
# A parallel readout to per-base error (phase 2) and DTW signal deviation
# (phase 2b): instead of averaging a signal over reads, each read is tested base by
# base against a canonical k-mer model built from the relationship's CONTROL, and
# the track is the fraction of reads called modified. Method: Rembo (P2-seq
# preprint), test from Tombo's model_sample_compare - see mod_calling_common.py.
#
# Two backends supply the per-read signal level; everything after that is shared:
#   uncalled4   dtw.current from the Uncalled4 BAM. Needs has_signal().
#   remora      Remora API on pod5 + the same BAM (it keeps the dorado move table,
#               so both backends score identical reads). Needs raw signal: the
#               sample's ^r pod5 path, or a '^pod5 <raw_sample> <path>' line for a
#               ^u sample. '^remora-levels <table>' turns on signal mapping refinement.
#
# The control is scored against a model built from itself, so reads are split in
# two by read name: the model uses half A, the control is scored on half B. The
# treatment is scored on all reads.
#
# Which control builds the model ({model} / {model_set}):
#   matched   the relationship's own control. Always built.
#   pooled    '^mod-pool <name> <raw_sample> ...' in CONFIG: one model from several
#             controls (replicates, other experiments). Built for every relationship
#             whose control is a member. Matched vs pooled reactivity is the
#             sample-to-sample reproducibility of the canonical signal.
#
# How a base is called ({thr}), both selectable at once for comparison:
#   p<P       --config mod_pval=P (default 0.02): Gaussian p-value below P.
#   fpr<F     --config mod_fpr=F: the cutoff a fraction F of the model's held-out
#             control reads exceed (rule perbase_mod_null measures it).
# Either key takes a comma-separated list.

import re
import sys

# A Snakefile generated before phase 2c existed defines none of these.
if "MOD_BACKENDS" not in globals():
    MOD_BACKENDS = []
if "POD5_PATHS" not in globals():
    POD5_PATHS = {}
if "REMORA_LEVELS" not in globals():
    REMORA_LEVELS = ""
if "MOD_POOLS" not in globals():
    MOD_POOLS = {}


def _threshold_list(key, default):
    value = config.get(key, default)
    return [] if value in (None, "") else [v.strip() for v in str(value).split(",")]


# Threshold tags used in output paths: p0.02, fpr0.02, ...
MOD_THRESHOLDS = ([f"p{v}" for v in _threshold_list("mod_pval", 0.02)]
                  + [f"fpr{v}" for v in _threshold_list("mod_fpr", None)])


def model_members(model):
    """Raw samples whose half-A reads build this model: a pool's members, or
    the one sample a matched model is named after."""
    return MOD_POOLS.get(model, [model])


def mod_pod5(raw_sample):
    """Raw signal for the remora backend: ^pod5 line first, else the ^r pod5 path."""
    return POD5_PATHS.get(raw_sample) or get_pod5_dir(raw_sample)


def has_mod_input(backend, raw_sample):
    """Return True if this backend can score this sample."""
    if not has_signal(raw_sample):
        return False
    return backend == "uncalled4" or mod_pod5(raw_sample) is not None


def mod_relationships(backend):
    """Relationship samples whose treatment and control both suit this backend."""
    return [s for s in SAMPLES
            if has_mod_input(backend, get_treatment(s))
            and has_mod_input(backend, get_control(s))]


def mod_model_sets(backend, sample):
    """'matched', plus every pool the relationship's control belongs to (and whose
    members this backend can all read)."""
    return ["matched"] + [
        name for name, members in MOD_POOLS.items()
        if get_control(sample) in members
        and all(has_mod_input(backend, m) for m in members)]


# A backend that was asked for and cannot run for a relationship must not just
# vanish from the DAG - that is how the signal branch was lost in js4022.
for _backend in MOD_BACKENDS:
    for _sample in SAMPLES:
        if _sample not in mod_relationships(_backend):
            print(f"WARNING: ^mod-calls {_backend} requested, but relationship "
                  f"'{_sample}' has no input for it and is SKIPPED "
                  f"(uncalled4 needs pod5 or ^u; remora also needs pod5: ^r or ^pod5).",
                  file=sys.stderr)
    for _name, _members in MOD_POOLS.items():
        _missing = [m for m in _members if not has_mod_input(_backend, m)]
        if _missing:
            print(f"WARNING: ^mod-pool {_name} is SKIPPED for backend {_backend}: "
                  f"no input for {', '.join(_missing)}.", file=sys.stderr)


def mod_sample_inputs(raw_sample, backend, genome):
    """BAM, index and (remora) pod5 of one raw sample."""
    supplied = get_uncalled4_bam(raw_sample)
    bam = supplied or f"data/uncalled4/{genome}/{raw_sample}.bam"
    pod5 = (mod_pod5(raw_sample) or []) if backend == "remora" else []
    return bam, bam + ".bai", pod5


def mod_inputs(wildcards):
    bam, bai, pod5 = mod_sample_inputs(wildcards.raw_sample, wildcards.backend,
                                       wildcards.genome)
    inputs = {"bam": bam, "bai": bai,
              "genome": f"resources/genomes/{wildcards.genome}.fa"}
    if wildcards.backend == "remora":
        inputs["pod5"] = pod5
        if REMORA_LEVELS:
            inputs["levels"] = REMORA_LEVELS
    return inputs


def mod_model_inputs(wildcards):
    """Every member of the model: one sample (matched) or several (^mod-pool)."""
    members = [mod_sample_inputs(m, wildcards.backend, wildcards.genome)
               for m in model_members(wildcards.model)]
    inputs = {"bam": [m[0] for m in members], "bai": [m[1] for m in members],
              "genome": f"resources/genomes/{wildcards.genome}.fa"}
    if wildcards.backend == "remora":
        inputs["pod5"] = [m[2] for m in members]
        if REMORA_LEVELS:
            inputs["levels"] = REMORA_LEVELS
    return inputs


def mod_backend_args(wildcards, input):
    """Backend-specific arguments shared by the model and the calling script."""
    if wildcards.backend != "remora":
        return ""
    args = f"--pod5 {input.pod5}"
    if REMORA_LEVELS:
        args += f" --levels {input.levels}"
    return args


def mod_half(wildcards):
    """A sample that helped build the model is scored on the reads it held out."""
    return "B" if wildcards.raw_sample in model_members(wildcards.model) else "all"


def mod_null_hists(wildcards):
    """Held-out histograms of every model member, both strands - the null for an
    fpr threshold. Nothing for a p-value threshold."""
    if not wildcards.thr.startswith("fpr"):
        return []
    return [f"data/perbase_mod_null/{wildcards.backend}/{wildcards.genome}/"
            f"{member}_vs_{wildcards.model}_{strand}.hist.tsv"
            for member in model_members(wildcards.model)
            for strand in STRANDS]


def mod_threshold_args(wildcards, input):
    if wildcards.thr.startswith("fpr"):
        return f"--fpr {wildcards.thr[3:]} --null-hist {input.null}"
    return f"--pval {wildcards.thr[1:]}"


MOD_SAMPLE_CONSTRAINT = "|".join(re.escape(rs) for rs in RAW_SAMPLES) or "NONE"
MOD_MODEL_CONSTRAINT = "|".join(
    re.escape(m) for m in list(RAW_SAMPLES) + list(MOD_POOLS)) or "NONE"
MOD_THR_CONSTRAINT = r"(p|fpr)[0-9.eE-]+"


rule kmer_signal_model:
    """Canonical per-k-mer signal model (n, mean, SD) from the half-A reads of one
    control, or of every member of a ^mod-pool."""
    input:
        unpack(mod_model_inputs)
    output:
        model="data/kmer_model/{backend}/{genome}/{model}.tsv"
    params:
        backend_args=mod_backend_args,
        mad=config.get("mod_mad", 15),
        max_obs=config.get("mod_max_obs_per_kmer", 20000)
    log:
        "logs/kmer_model/{backend}/{genome}/{model}.log"
    benchmark:
        "benchmarks/phase2c/kmer_signal_model/{backend}/{genome}/{model}.tsv"
    wildcard_constraints:
        backend="uncalled4|remora",
        model=MOD_MODEL_CONSTRAINT
    shell:
        """
        python3 workflow/scripts/kmer_signal_model.py \
            --backend {wildcards.backend} \
            -i {input.bam} \
            -g {input.genome} \
            {params.backend_args} \
            --half A \
            --mad {params.mad} \
            --max-obs-per-kmer {params.max_obs} \
            -o {output.model} \
            2>&1 | tee {log}
        """


rule perbase_mod_null:
    """Null distribution of the calling statistic: a model member's held-out reads
    scored against the model. Only requested for fpr thresholds."""
    input:
        unpack(mod_inputs),
        model="data/kmer_model/{backend}/{genome}/{model}.tsv"
    output:
        hist="data/perbase_mod_null/{backend}/{genome}/{raw_sample}_vs_{model}_{strand}.hist.tsv"
    params:
        backend_args=mod_backend_args,
        min_obs=config.get("mod_min_kmer_obs", 30),
        fisher_lag=config.get("mod_fisher_lag", 0)
    log:
        "logs/perbase_mod_null/{backend}/{genome}/{raw_sample}_vs_{model}_{strand}.log"
    benchmark:
        "benchmarks/phase2c/perbase_mod_null/{backend}/{genome}/{raw_sample}_vs_{model}_{strand}.tsv"
    wildcard_constraints:
        backend="uncalled4|remora",
        strand="for|rev",
        raw_sample=MOD_SAMPLE_CONSTRAINT,
        model=MOD_MODEL_CONSTRAINT
    shell:
        """
        python3 workflow/scripts/perbase_mod_calls.py \
            --backend {wildcards.backend} \
            -i {input.bam} \
            -g {input.genome} \
            {params.backend_args} \
            -m {input.model} \
            -s {wildcards.strand} \
            --half B \
            --min-kmer-obs {params.min_obs} \
            --fisher-lag {params.fisher_lag} \
            --hist-out {output.hist} \
            2>&1 | tee {log}
        """


rule perbase_mod_calls:
    """Score one sample strand against a k-mer model; modified fraction per position.
    A sample that is part of the model is scored on its held-out half B.
    """
    input:
        unpack(mod_inputs),
        model="data/kmer_model/{backend}/{genome}/{model}.tsv",
        null=mod_null_hists
    output:
        calls="data/perbase_mod/{backend}/{genome}/{thr}/{raw_sample}_vs_{model}_{strand}.txt.gz"
    params:
        backend_args=mod_backend_args,
        threshold_args=mod_threshold_args,
        half=mod_half,
        min_obs=config.get("mod_min_kmer_obs", 30),
        fisher_lag=config.get("mod_fisher_lag", 0),
        # Per-read table (read_id, pos, z, p) for single-read heatmaps. Not a tracked
        # output: it is optional and can be very large on anything but model oligos.
        # z and p do not depend on the threshold, so there is no {thr} in its path.
        per_read=lambda w: (
            f"--per-read data/perbase_mod_reads/{w.backend}/{w.genome}/"
            f"{w.raw_sample}_vs_{w.model}_{w.strand}.tsv.gz"
            if config.get("mod_per_read", False) else "")
    log:
        "logs/perbase_mod/{backend}/{genome}/{thr}/{raw_sample}_vs_{model}_{strand}.log"
    benchmark:
        "benchmarks/phase2c/perbase_mod_calls/{backend}/{genome}/{thr}/{raw_sample}_vs_{model}_{strand}.tsv"
    wildcard_constraints:
        backend="uncalled4|remora",
        strand="for|rev",
        thr=MOD_THR_CONSTRAINT,
        raw_sample=MOD_SAMPLE_CONSTRAINT,
        model=MOD_MODEL_CONSTRAINT
    shell:
        """
        python3 workflow/scripts/perbase_mod_calls.py \
            --backend {wildcards.backend} \
            -i {input.bam} \
            -g {input.genome} \
            {params.backend_args} \
            -m {input.model} \
            -s {wildcards.strand} \
            --half {params.half} \
            {params.threshold_args} \
            --min-kmer-obs {params.min_obs} \
            --fisher-lag {params.fisher_lag} \
            {params.per_read} \
            -o {output.calls} \
            2>&1 | tee {log}
        """
