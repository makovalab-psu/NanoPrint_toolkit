# Modification-calling parameter sweep
#
# Re-runs direct modification calling (phases 2c/3c) over a table of parameter
# sets - alignment filter, model options, test options, threshold - and tabulates
# every set side by side. Meant for small datasets (model oligos, plasmids) where
# the defaults can be tuned cheaply before they are applied to a genome.
#
# It sits on top of a finished main workflow and does not touch its outputs:
# everything is written under sweep/. It reads the main workflow's
# data/aligned_reads/ BAMs, pod5 paths, relationships and ^mod-pool definitions.
#
#   ./workflow/scripts/CONFIG.sh -i CONFIG -o Snakefile      # as usual
#   snakemake -s workflow/sweeps/mod_sweep.smk --cores 30 --keep-going \
#       --rerun-triggers mtime \
#       --config sweep_table=workflow/sweeps/mod_sweep_example.tsv
#
# --rerun-triggers mtime: this file includes the main workflow, so without it a
# main rule whose code or params changed since its outputs were built (as
# filter_alignments did when filter_mapq was added) is rebuilt along with everything
# downstream of it, although the sweep only wants to read those files.
#
# Sweep table: tab-separated, one parameter set per row, '#' lines ignored.
#   set           name of the parameter set (unique; used in the output tables)
#   backend       uncalled4 | remora
#   mapq          minimum MAPQ of the alignment filter; 0 = every mapped primary read
#   mad           MAD multiple of the model's outlier filter
#   max_obs       observations kept per k-mer in the model (0 = all)
#   flank         error-mask flank of the model: default (= k-1) | N | off
#   scale         sd (k-mer mean and SD) | mad (median and 1.4826 x MAD)
#   min_kmer_obs  k-mers the model saw fewer times than this are not scored
#   fisher_lag    Fisher's-method window along the read (0 = off)
#   threshold     p<value> or fpr<value>, e.g. p0.02, fpr0.02
#
# Optional --config keys:
#   sweep_name=<name>             output folder under sweep/tables/ (default: table name)
#   sweep_relationships=a,b       only these relationship samples (default: all)
#   sweep_genomes=g1,g2           only these genomes (default: all)
#   sweep_pools=true              also score against every ^mod-pool the control is in
#   filter_mapq=<N>               the MAPQ the MAIN workflow was run with, if it was
#                                 given by --config and not by ^filter-mapq;
#                                 sets of that MAPQ reuse its Uncalled4 BAMs
#   sweep_min_cov=<N>             coverage both samples need for a position to enter
#                                 the per-position summary columns (default:
#                                 reactivity_cov_threshold, 10). The pooled columns
#                                 use every position both samples scored.
#   sweep_position_rows=auto|all|none
#                                 rows of positions.tsv.gz (default auto: all for
#                                 model oligos, none once a call file has more than
#                                 200,000 positions - a genome)
#   sweep_regions=<bed>           positions.tsv.gz only inside these intervals, and
#                                 an extra scope=regions summary row per contig.
#                                 With sweep_position_rows=none: the summary rows
#                                 only, for sites too thinly covered to read singly
#   sweep_table_threads=<N>       processes for the table (default 8; ~1.5 GB each
#                                 on a 4.2 Mb genome)
#   sweep_read_summary_sets=a,b   also score these sets one row per read
#                                 (perbase_mod_calls.py --read-summary) and compare
#                                 treated with control per read: are the calls
#                                 concentrated in some molecules? Off by default.
#   sweep_read_min_scored=<N>     reads with fewer scored positions are left out of
#                                 that comparison (default 1000)
#
# Genomes. A ^u sample (already Uncalled4-aligned) has no pre-filter reads, so a set
# BELOW the main workflow's MAPQ cannot include it. A set ABOVE it can: the supplied
# BAM is filtered to the higher MAPQ (sweep_subset_uncalled4) - the DTW of a read
# does not depend on which other reads were aligned with it. The model cap (max_obs)
# matters on a genome: all 262,144 9-mers are present, so use ~1000, not the 20,000
# chosen on oligos (js4033).
#
# Outputs:
#   sweep/tables/<name>/positions.tsv.gz   one row per set x relationship x position
#                                          (see sweep_position_rows, sweep_regions)
#   sweep/tables/<name>/summary.tsv        one row per set x relationship x contig
#                                          (x scope, with sweep_regions)
#   sweep/tables/<name>/manifest.tsv       which files each row came from
#   sweep/tables/<name>/reads_summary.tsv  with sweep_read_summary_sets: per-read
#   sweep/tables/<name>/reads_hist.tsv     call rates, treated vs control
#   sweep/tables/<name>/reads_manifest.tsv (read_summary_table.py)
#
# Intermediate files are keyed by their parameters, so sets that share a filter or
# a model share the files, and a second sweep reuses whatever the first one built.

import os
import csv
import re

SWEEP_TABLE = config.get("sweep_table", "workflow/sweeps/mod_sweep_example.tsv")
SWEEP_NAME = config.get(
    "sweep_name", os.path.splitext(os.path.basename(SWEEP_TABLE))[0])
SWEEP_COLUMNS = ["set", "backend", "mapq", "mad", "max_obs", "flank", "scale",
                 "min_kmer_obs", "fisher_lag", "threshold"]


def read_sweep_table(path):
    with open(path) as fh:
        lines = [line for line in fh if line.strip() and not line.startswith("#")]
    rows = list(csv.DictReader(lines, delimiter="\t"))
    if not rows:
        raise ValueError(f"Sweep table {path} has no parameter sets.")
    missing = [c for c in SWEEP_COLUMNS if c not in rows[0]]
    if missing:
        raise ValueError(f"Sweep table {path} is missing column(s): {', '.join(missing)}. "
                         f"Columns are tab-separated: {' '.join(SWEEP_COLUMNS)}")
    seen = set()
    for row in rows:
        for key in SWEEP_COLUMNS:
            row[key] = (row[key] or "").strip()
            if not row[key]:
                raise ValueError(f"Sweep table {path}: set '{row['set']}' has no {key}.")
        if row["set"] in seen:
            raise ValueError(f"Sweep table {path}: set name '{row['set']}' is used twice.")
        seen.add(row["set"])
        if row["backend"] not in ("uncalled4", "remora"):
            raise ValueError(f"Set '{row['set']}': backend must be uncalled4 or remora.")
        if row["scale"] not in ("sd", "mad"):
            raise ValueError(f"Set '{row['set']}': scale must be sd or mad.")
        if not re.fullmatch(r"default|off|\d+", row["flank"]):
            raise ValueError(f"Set '{row['set']}': flank must be default, off or a number.")
        if not re.fullmatch(r"(p|fpr)[0-9.eE-]+", row["threshold"]):
            raise ValueError(f"Set '{row['set']}': threshold must look like p0.02 or fpr0.02.")
    return rows


SWEEP_SETS = read_sweep_table(SWEEP_TABLE)


def _config_list(key):
    value = config.get(key)
    return [v.strip() for v in str(value).split(",")] if value not in (None, "") else None


def sweep_model_tag(row):
    return f"mad{row['mad']}_cap{row['max_obs']}_flank{row['flank']}"


def sweep_test_tag(row):
    return f"min{row['min_kmer_obs']}_fl{row['fisher_lag']}_{row['scale']}"


def sweep_units():
    """Every (parameter set, genome, relationship, control model, strand) to score."""
    relationships = _config_list("sweep_relationships") or SAMPLES
    unknown = [r for r in relationships if r not in SAMPLES]
    if unknown:
        raise ValueError(f"sweep_relationships: not in CONFIG: {', '.join(unknown)}")
    genomes = _config_list("sweep_genomes") or GENOMES
    use_pools = str(config.get("sweep_pools", "false")).lower() in ("true", "1", "yes")
    units = []
    for row in SWEEP_SETS:
        for genome in genomes:
            for sample in relationships:
                treatment, control = get_treatment(sample), get_control(sample)
                if not (has_mod_input(row["backend"], treatment)
                        and has_mod_input(row["backend"], control)):
                    continue
                model_sets = mod_model_sets(row["backend"], sample) if use_pools else ["matched"]
                for model_set in model_sets:
                    model = control if model_set == "matched" else model_set
                    base = (f"sweep/perbase_mod/{row['backend']}/mapq{row['mapq']}/"
                            f"{sweep_model_tag(row)}/{row['threshold']}_{sweep_test_tag(row)}/"
                            f"{genome}")
                    for strand in STRANDS:
                        units.append(dict(
                            row, genome=genome, relationship=sample, model_set=model_set,
                            model=model, strand=strand,
                            treatment_file=f"{base}/{treatment}_vs_{model}_{strand}.txt.gz",
                            control_file=f"{base}/{control}_vs_{model}_{strand}.txt.gz"))
    if not units:
        raise ValueError("The sweep selects nothing to score. Check sweep_relationships, "
                         "sweep_genomes, and that the samples have signal data (and pod5 "
                         "for remora).")
    return units


def sweep_call_files(wildcards):
    files = []
    for unit in sweep_units():
        files += [unit["treatment_file"], unit["control_file"]]
    return sorted(set(files))


SWEEP_READ_SETS = _config_list("sweep_read_summary_sets") or []
_unknown_read_sets = sorted(set(SWEEP_READ_SETS) - {row["set"] for row in SWEEP_SETS})
if _unknown_read_sets:
    raise ValueError(f"sweep_read_summary_sets: not in {SWEEP_TABLE}: "
                     f"{', '.join(_unknown_read_sets)}")


def read_summary_file(call_file):
    """The --read-summary file scored alongside a perbase_mod call file."""
    return ("sweep/read_summary/" + call_file[len("sweep/perbase_mod/"):]
            ).replace(".txt.gz", ".reads.tsv.gz")


def sweep_read_units():
    return [unit for unit in sweep_units() if unit["set"] in SWEEP_READ_SETS]


def sweep_read_files(wildcards):
    files = []
    for unit in sweep_read_units():
        files += [read_summary_file(unit["treatment_file"]),
                  read_summary_file(unit["control_file"])]
    return sorted(set(files))


# Defined before the main Snakefile is included so that it is the default target.
rule sweep_all:
    input:
        f"sweep/tables/{SWEEP_NAME}/summary.tsv",
        f"sweep/tables/{SWEEP_NAME}/positions.tsv.gz",
        [f"sweep/tables/{SWEEP_NAME}/reads_summary.tsv"] if SWEEP_READ_SETS else []


# The main workflow: sample names, relationships, pod5 paths, pools, and the
# helper functions of phases 0 and 2c. Its 'rule all' is not the target here.
include: os.path.abspath(config.get("main_snakefile", "Snakefile"))

if REMORA_LEVELS == "" and any(row["backend"] == "remora" for row in SWEEP_SETS):
    raise ValueError("The sweep table has remora sets but CONFIG has no ^remora-levels "
                     "line. The remora backend is not usable without the level table.")

# MAPQ the main workflow filtered at. Sets at this MAPQ reuse its Uncalled4 BAMs
# instead of re-aligning signal, which also makes them the exact baseline.
MAIN_FILTER_MAPQ = int(config.get("filter_mapq", FILTER_MAPQ))


def sweep_sample_bam(raw_sample, genome, mapq):
    """Uncalled4 BAM of one sample at one filter setting."""
    supplied = get_uncalled4_bam(raw_sample)
    if int(mapq) == MAIN_FILTER_MAPQ:
        return supplied or f"data/uncalled4/{genome}/{raw_sample}.bam"
    if supplied:
        if int(mapq) < MAIN_FILTER_MAPQ:
            raise ValueError(
                f"{raw_sample} was supplied already Uncalled4-aligned (^u) at MAPQ >= "
                f"{MAIN_FILTER_MAPQ}, so its lower-MAPQ reads do not exist and it cannot "
                f"be in a mapq={mapq} set. Sets at mapq >= {MAIN_FILTER_MAPQ} can include it.")
        return f"sweep/uncalled4_subset/mapq{mapq}/{genome}/{raw_sample}.bam"
    return f"sweep/uncalled4/mapq{mapq}/{genome}/{raw_sample}.bam"


def sweep_inputs(raw_samples, wildcards):
    bams = [sweep_sample_bam(rs, wildcards.genome, wildcards.mapq) for rs in raw_samples]
    inputs = {"bam": bams, "bai": [b + ".bai" for b in bams],
              "genome": f"resources/genomes/{wildcards.genome}.fa"}
    if wildcards.backend == "remora":
        inputs["pod5"] = [mod_pod5(rs) or [] for rs in raw_samples]
        inputs["levels"] = REMORA_LEVELS
    return inputs


def sweep_flank_arg(wildcards):
    if wildcards.flank == "default":
        return ""
    return "--clean-flank -1" if wildcards.flank == "off" else f"--clean-flank {wildcards.flank}"


SWEEP_MODEL_DIR = "sweep/kmer_model/{backend}/mapq{mapq}/mad{mad}_cap{cap}_flank{flank}/{genome}"
SWEEP_CALL_DIR = ("{backend}/mapq{mapq}/mad{mad}_cap{cap}_flank{flank}/"
                  "{thr}_min{min_obs}_fl{lag}_{scale}/{genome}")
SWEEP_NULL_DIR = ("{backend}/mapq{mapq}/mad{mad}_cap{cap}_flank{flank}/"
                  "min{min_obs}_fl{lag}_{scale}/{genome}")

wildcard_constraints:
    mapq=r"\d+",
    mad=r"[0-9.]+",
    cap=r"\d+",
    flank=r"default|off|\d+",
    min_obs=r"\d+",
    lag=r"\d+",
    scale="sd|mad"


rule sweep_filter_alignments:
    """Re-filter the main workflow's aligned reads at another minimum MAPQ."""
    input:
        bam="data/aligned_reads/{genome}/{raw_sample}.bam"
    output:
        bam="sweep/filtered_alignments/mapq{mapq}/{genome}/{raw_sample}.bam"
    log:
        "logs/sweep/filter_alignments/mapq{mapq}/{genome}/{raw_sample}.log"
    benchmark:
        "benchmarks/sweep/filter_alignments/mapq{mapq}/{genome}/{raw_sample}.tsv"
    wildcard_constraints:
        raw_sample=MOD_SAMPLE_CONSTRAINT
    shell:
        """
        workflow/scripts/Filter_alignments.sh \
            -i {input.bam} \
            -o {output.bam} \
            -q {wildcards.mapq} \
            2>&1 | tee {log}
        """


rule sweep_uncalled4_align:
    """Signal-align the re-filtered reads (same script and settings as phase 0)."""
    input:
        bam="sweep/filtered_alignments/mapq{mapq}/{genome}/{raw_sample}.bam",
        pod5=lambda wildcards: pod5_input(wildcards.raw_sample),
        genome="resources/genomes/{genome}.fa"
    output:
        bam="sweep/uncalled4/mapq{mapq}/{genome}/{raw_sample}.bam",
        bai="sweep/uncalled4/mapq{mapq}/{genome}/{raw_sample}.bam.bai"
    log:
        "logs/sweep/uncalled4_align/mapq{mapq}/{genome}/{raw_sample}.log"
    benchmark:
        "benchmarks/sweep/uncalled4_align/mapq{mapq}/{genome}/{raw_sample}.tsv"
    wildcard_constraints:
        raw_sample=MOD_SAMPLE_CONSTRAINT
    threads: 8
    shell:
        """
        workflow/scripts/Uncalled4_align.sh \
            -i {input.bam} \
            -p {input.pod5} \
            -g {input.genome} \
            -o {output.bam} \
            -t {threads} \
            2>&1 | tee {log}
        """


rule sweep_subset_uncalled4:
    """Raise the MAPQ floor of a supplied (^u) Uncalled4 BAM; DTW tags are kept."""
    input:
        bam=lambda wildcards: get_uncalled4_bam(wildcards.raw_sample),
        bai=lambda wildcards: get_uncalled4_bam(wildcards.raw_sample) + ".bai"
    output:
        bam="sweep/uncalled4_subset/mapq{mapq}/{genome}/{raw_sample}.bam",
        bai="sweep/uncalled4_subset/mapq{mapq}/{genome}/{raw_sample}.bam.bai"
    log:
        "logs/sweep/subset_uncalled4/mapq{mapq}/{genome}/{raw_sample}.log"
    benchmark:
        "benchmarks/sweep/subset_uncalled4/mapq{mapq}/{genome}/{raw_sample}.tsv"
    wildcard_constraints:
        raw_sample="|".join(re.escape(rs) for rs in UNCALLED4_PATHS) or "NONE"
    threads: 4
    shell:
        """
        workflow/scripts/Subset_uncalled4_mapq.sh \
            -i {input.bam} \
            -o {output.bam} \
            -q {wildcards.mapq} \
            -t {threads} \
            2>&1 | tee {log}
        """


rule sweep_kmer_signal_model:
    """K-mer model for one filter setting and one set of model options."""
    input:
        unpack(lambda wildcards: sweep_inputs(model_members(wildcards.model), wildcards))
    output:
        model=SWEEP_MODEL_DIR + "/{model}.tsv"
    params:
        backend_args=mod_backend_args,
        flank=sweep_flank_arg
    log:
        "logs/" + SWEEP_MODEL_DIR + "/{model}.log"
    benchmark:
        "benchmarks/" + SWEEP_MODEL_DIR + "/{model}.tsv"
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
            --mad {wildcards.mad} \
            --max-obs-per-kmer {wildcards.cap} \
            {params.flank} \
            -o {output.model} \
            2>&1 | tee {log}
        """


rule sweep_perbase_mod_null:
    """Held-out null of the calling statistic, for fpr thresholds."""
    input:
        unpack(lambda wildcards: sweep_inputs([wildcards.raw_sample], wildcards)),
        model=SWEEP_MODEL_DIR + "/{model}.tsv"
    output:
        hist="sweep/perbase_mod_null/" + SWEEP_NULL_DIR + "/{raw_sample}_vs_{model}_{strand}.hist.tsv"
    params:
        backend_args=mod_backend_args
    log:
        "logs/sweep/perbase_mod_null/" + SWEEP_NULL_DIR + "/{raw_sample}_vs_{model}_{strand}.log"
    benchmark:
        "benchmarks/sweep/perbase_mod_null/" + SWEEP_NULL_DIR + "/{raw_sample}_vs_{model}_{strand}.tsv"
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
            --scale {wildcards.scale} \
            --min-kmer-obs {wildcards.min_obs} \
            --fisher-lag {wildcards.lag} \
            --hist-out {output.hist} \
            2>&1 | tee {log}
        """


def sweep_null_hists(wildcards):
    if not wildcards.thr.startswith("fpr"):
        return []
    return expand(
        "sweep/perbase_mod_null/" + SWEEP_NULL_DIR + "/{member}_vs_{model}_{strand}.hist.tsv",
        member=model_members(wildcards.model), strand=STRANDS,
        **{k: getattr(wildcards, k) for k in ("backend", "mapq", "mad", "cap", "flank",
                                              "min_obs", "lag", "scale", "genome",
                                              "model")})


rule sweep_perbase_mod_calls:
    """Score one sample strand under one parameter set."""
    input:
        unpack(lambda wildcards: sweep_inputs([wildcards.raw_sample], wildcards)),
        model=SWEEP_MODEL_DIR + "/{model}.tsv",
        null=sweep_null_hists
    output:
        calls="sweep/perbase_mod/" + SWEEP_CALL_DIR + "/{raw_sample}_vs_{model}_{strand}.txt.gz"
    params:
        backend_args=mod_backend_args,
        threshold_args=mod_threshold_args,
        half=mod_half
    log:
        "logs/sweep/perbase_mod/" + SWEEP_CALL_DIR + "/{raw_sample}_vs_{model}_{strand}.log"
    benchmark:
        "benchmarks/sweep/perbase_mod/" + SWEEP_CALL_DIR + "/{raw_sample}_vs_{model}_{strand}.tsv"
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
            --scale {wildcards.scale} \
            --min-kmer-obs {wildcards.min_obs} \
            --fisher-lag {wildcards.lag} \
            -o {output.calls} \
            2>&1 | tee {log}
        """


rule sweep_table:
    """Tabulate every parameter set: per position, and summarised per contig."""
    input:
        calls=sweep_call_files,
        table=SWEEP_TABLE,
        regions=config.get("sweep_regions") or []
    output:
        manifest="sweep/tables/{name}/manifest.tsv",
        positions="sweep/tables/{name}/positions.tsv.gz",
        summary="sweep/tables/{name}/summary.tsv"
    params:
        cov=config.get("sweep_min_cov", config.get("reactivity_cov_threshold", 10)),
        rows=config.get("sweep_position_rows", "auto"),
        regions=lambda wildcards, input: f"--regions {input.regions}" if input.regions else ""
    log:
        "logs/sweep/tables/{name}.log"
    threads: int(config.get("sweep_table_threads", 8))
    run:
        units = sweep_units()
        columns = SWEEP_COLUMNS + ["genome", "relationship", "model_set", "model",
                                   "strand", "treatment_file", "control_file"]
        with open(output.manifest, "w") as fh:
            fh.write("\t".join(columns) + "\n")
            for unit in units:
                fh.write("\t".join(str(unit[c]) for c in columns) + "\n")
        shell(
            "python3 workflow/scripts/mod_sweep_table.py "
            "-m {output.manifest} -p {output.positions} -s {output.summary} "
            "-c {params.cov} -t {threads} --position-rows {params.rows} {params.regions} "
            "2>&1 | tee {log}")


# Per-read call rates (sweep_read_summary_sets). A separate scoring pass, not an
# extra output of sweep_perbase_mod_calls: adding an output there would rerun every
# set of a finished sweep. Same inputs and options, so the calls are the same.
rule sweep_read_summary:
    """One row per read of one sample strand under one parameter set."""
    input:
        unpack(lambda wildcards: sweep_inputs([wildcards.raw_sample], wildcards)),
        model=SWEEP_MODEL_DIR + "/{model}.tsv",
        null=sweep_null_hists
    output:
        reads="sweep/read_summary/" + SWEEP_CALL_DIR + "/{raw_sample}_vs_{model}_{strand}.reads.tsv.gz"
    params:
        backend_args=mod_backend_args,
        threshold_args=mod_threshold_args,
        half=mod_half
    log:
        "logs/sweep/read_summary/" + SWEEP_CALL_DIR + "/{raw_sample}_vs_{model}_{strand}.log"
    benchmark:
        "benchmarks/sweep/read_summary/" + SWEEP_CALL_DIR + "/{raw_sample}_vs_{model}_{strand}.tsv"
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
            --scale {wildcards.scale} \
            --min-kmer-obs {wildcards.min_obs} \
            --fisher-lag {wildcards.lag} \
            --read-summary {output.reads} \
            2>&1 | tee {log}
        """


rule sweep_read_table:
    """Per-read call rates of every sweep_read_summary_sets set, treated vs control."""
    input:
        reads=sweep_read_files,
        table=SWEEP_TABLE
    output:
        manifest="sweep/tables/{name}/reads_manifest.tsv",
        summary="sweep/tables/{name}/reads_summary.tsv",
        hist="sweep/tables/{name}/reads_hist.tsv"
    params:
        min_scored=config.get("sweep_read_min_scored", 1000)
    log:
        "logs/sweep/tables/{name}_reads.log"
    run:
        columns = SWEEP_COLUMNS + ["genome", "relationship", "model_set", "model", "strand"]
        with open(output.manifest, "w") as fh:
            fh.write("\t".join(columns + ["treatment_reads", "control_reads"]) + "\n")
            for unit in sweep_read_units():
                fh.write("\t".join([str(unit[c]) for c in columns]
                                   + [read_summary_file(unit["treatment_file"]),
                                      read_summary_file(unit["control_file"])]) + "\n")
        shell(
            "python3 workflow/scripts/read_summary_table.py "
            "-m {output.manifest} -s {output.summary} --hist {output.hist} "
            "--min-scored {params.min_scored} "
            "2>&1 | tee {log}")
