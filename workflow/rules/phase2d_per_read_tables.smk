# Phase 2d: Per-Read, Per-Nucleotide Tables (optional)
# Off by default. Turned on in CONFIG, one line per table wanted:
#
#   ^per-read	uncalled4     signal of every read at every position, from the Uncalled4 BAM
#   ^per-read	remora        the same through the Remora API (needs ^remora-levels, pod5)
#   ^per-read	mod-calls     z, p and the modified/unmodified call of every read at every
#                             position, for each ^mod-calls backend
#
# These are the single-molecule view: nothing is averaged, and every table carries
# the read name, so one molecule can be followed across all three.
#
# They are one row per aligned base per read - practical for model oligos, plasmids
# and single loci, not for a genome - which is why nothing here is built unless
# CONFIG asks for it.
#
#   data/perread_signal/{backend}/{genome}/{raw_sample}_{strand}.tsv.gz
#   data/perread_mod/{backend}/{genome}/{thr}/{raw_sample}_vs_{model}_{strand}.tsv.gz
#
# The uncalled4 table is written by perread_signal.py from the BAM's DTW tags, not
# by 'uncalled4 convert': the convert TSV names a read only by a running number and
# cannot be joined to the other two. uncalled4_convert_tsv (phase 0) is unchanged and
# can still be requested by name.

# A Snakefile generated before phase 2d existed does not define this.
if "PER_READ" not in globals():
    PER_READ = []


def per_read_targets():
    """Everything the ^per-read lines ask for."""
    targets = []
    for backend in ("uncalled4", "remora"):
        if backend in PER_READ:
            targets += [
                f"data/perread_signal/{backend}/{genome}/{rs}_{strand}.tsv.gz"
                for genome in GENOMES
                for rs in RAW_SAMPLES if has_mod_input(backend, rs)
                for strand in STRANDS]
    if "mod-calls" in PER_READ:
        for backend in MOD_BACKENDS:
            for sample in mod_relationships(backend):
                treatment, control = get_treatment(sample), get_control(sample)
                for model_set in mod_model_sets(backend, sample):
                    model = control if model_set == "matched" else model_set
                    targets += [
                        f"data/perread_mod/{backend}/{genome}/{thr}/"
                        f"{rs}_vs_{model}_{strand}.tsv.gz"
                        for genome in GENOMES
                        for thr in MOD_THRESHOLDS
                        for rs in (treatment, control)
                        for strand in STRANDS]
    return sorted(set(targets))


# Same reasoning as the ^mod-calls warning in phase 2c: a table that was asked for
# and cannot be built for a sample must say so rather than be missing.
for _backend in ("uncalled4", "remora"):
    if _backend in PER_READ:
        for _rs in RAW_SAMPLES:
            if not has_mod_input(_backend, _rs):
                print(f"WARNING: ^per-read {_backend} requested, but sample '{_rs}' has "
                      f"no input for it and is SKIPPED (uncalled4 needs pod5 or ^u; "
                      f"remora also needs pod5: ^r or ^pod5).", file=sys.stderr)


rule perread_signal:
    """Per-read, per-nucleotide signal table for one sample strand."""
    input:
        unpack(mod_inputs)
    output:
        tsv="data/perread_signal/{backend}/{genome}/{raw_sample}_{strand}.tsv.gz"
    params:
        backend_args=mod_backend_args
    log:
        "logs/perread_signal/{backend}/{genome}/{raw_sample}_{strand}.log"
    benchmark:
        "benchmarks/phase2d/perread_signal/{backend}/{genome}/{raw_sample}_{strand}.tsv"
    wildcard_constraints:
        backend="uncalled4|remora",
        strand="for|rev",
        raw_sample=MOD_SAMPLE_CONSTRAINT
    shell:
        """
        python3 workflow/scripts/perread_signal.py \
            --backend {wildcards.backend} \
            -i {input.bam} \
            -g {input.genome} \
            {params.backend_args} \
            -s {wildcards.strand} \
            -o {output.tsv} \
            2>&1 | tee {log}
        """


rule perread_mod_calls:
    """Per-read modification calls for one sample strand: z, p and the call of every
    read at every scored position. Same model, reads, threshold and settings as
    perbase_mod_calls, whose per-position table is the average of this one.
    """
    input:
        unpack(mod_inputs),
        model="data/kmer_model/{backend}/{genome}/{model}.tsv",
        null=mod_null_hists
    output:
        tsv="data/perread_mod/{backend}/{genome}/{thr}/{raw_sample}_vs_{model}_{strand}.tsv.gz"
    params:
        backend_args=mod_backend_args,
        threshold_args=mod_threshold_args,
        half=mod_half,
        min_obs=config.get("mod_min_kmer_obs", 30),
        fisher_lag=MOD_FISHER_LAG
    log:
        "logs/perread_mod/{backend}/{genome}/{thr}/{raw_sample}_vs_{model}_{strand}.log"
    benchmark:
        "benchmarks/phase2d/perread_mod_calls/{backend}/{genome}/{thr}/{raw_sample}_vs_{model}_{strand}.tsv"
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
            --per-read {output.tsv} \
            2>&1 | tee {log}
        """
