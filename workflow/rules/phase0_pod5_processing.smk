# Phase 0: Pod5 Raw Data Processing
# Dorado basecalling and Uncalled4 signal alignment.
# Rules here only execute when pod5 input is detected for a sample.
#
# A sample enters the workflow in one of three ways, decided by its CONFIG line:
#   ^r + pod5 path   basecall -> map -> filter -> uncalled4 align (phase 0 + 1)
#   ^r + reads/BAM   map -> filter (phase 1); no signal data, so no phase 2b/3b/4b
#   ^u + BAM         already Uncalled4-aligned; phases 0 and 1 do not run at all
#
# has_pod5(), has_uncalled4(), has_signal() and the path helpers are defined here
# so they are available to all later phases (this file is included first).

import os


# A Snakefile generated before ^u existed does not define UNCALLED4_PATHS. Such a
# workflow has no ^u samples by definition, so an empty mapping is the right answer
# and the old Snakefile keeps working instead of dying on a NameError several
# includes later.
if "UNCALLED4_PATHS" not in globals():
    UNCALLED4_PATHS = {}


def get_pod5_dir(raw_sample):
    """Return the pod5 source path for a sample, or None if not pod5 input.

    The path comes from RAW_PATHS (set by CONFIG.sh from the ^r line).
    A single .pod5 file is returned as-is.
    A directory is returned as-is if any .pod5 file exists anywhere under it
    (os.walk recurses into subdirectories so nested sequencer output trees work).
    """
    path = RAW_PATHS.get(raw_sample)
    if path is None:
        return None
    if os.path.isfile(path) and path.endswith(".pod5"):
        return path
    if os.path.isdir(path):
        for _root, _dirs, files in os.walk(path):
            if any(f.endswith(".pod5") for f in files):
                return path  # top-level dir; Dorado_basecall.sh adds --recursive, Uncalled4_align.sh uses find
    return None


def has_pod5(raw_sample):
    """Return True if this sample's raw input is pod5 data."""
    return get_pod5_dir(raw_sample) is not None


def pod5_input(raw_sample):
    """Pod5 input for a rule, or no input at all when the sample is not pod5.

    Snakemake instantiates the job that produces a file even when that file
    already exists, purely to decide whether it is up to date — and doing so
    evaluates the input functions. Returning None there raises
    "Function did not return str or iterable of str. Encountered: [None]"
    during DAG construction, before the up-to-date check can happen. That
    breaks any workflow where data/uncalled4/ or data/basecalled/ was supplied
    from outside the pipeline rather than produced by it.

    An empty list means "this rule has no pod5 input", which is correct for a
    non-pod5 sample: the rule is not applicable and must never run.
    """
    return get_pod5_dir(raw_sample) or []


def get_uncalled4_bam(raw_sample):
    """Return the supplied Uncalled4 BAM for a sample (^u in CONFIG), or None.

    UNCALLED4_PATHS is written by CONFIG.sh, which has already checked that the
    path ends in .bam, that the sample is not also declared by ^r, and that an
    index exists next to it.
    """
    return UNCALLED4_PATHS.get(raw_sample)


def has_uncalled4(raw_sample):
    """Return True if this sample was supplied already Uncalled4-aligned."""
    return get_uncalled4_bam(raw_sample) is not None


def has_signal(raw_sample):
    """Return True if DTW signal data is or will be available for this sample.

    Either there is pod5 to align (^r) or an Uncalled4 BAM was supplied (^u).
    This is what phases 2b/3b/4b need; they read DTW tags, not raw signal, so a
    supplied BAM serves them exactly as a freshly aligned one does.
    """
    return has_pod5(raw_sample) or has_uncalled4(raw_sample)


def uncalled4_bam(wildcards):
    """The Uncalled4 BAM for a sample: supplied by ^u, or built by uncalled4_align.

    Returning the supplied path verbatim is what makes the DAG start after signal
    alignment — nothing upstream of it is ever requested, so no rule that would
    rebuild it is instantiated.
    """
    supplied = get_uncalled4_bam(wildcards.raw_sample)
    if supplied:
        return supplied
    return f"data/uncalled4/{wildcards.genome}/{wildcards.raw_sample}.bam"


def uncalled4_bai(wildcards):
    """Index of the Uncalled4 BAM; alongside the BAM in both cases."""
    return uncalled4_bam(wildcards) + ".bai"


def igv_sources(raw_sample):
    """IGV export sources applicable to a sample.

    A ^u sample has no pre-filter BAM — those reads were dropped before Uncalled4
    ran — so only the filtered_alignments export exists for it.
    """
    if has_uncalled4(raw_sample):
        return ["filtered_alignments"]
    return IGV_SOURCES


rule dorado_basecall:
    """Basecall pod5 files with Dorado (emits move tables for Uncalled4)."""
    input:
        pod5=lambda wildcards: pod5_input(wildcards.raw_sample)
    output:
        bam="data/basecalled/{raw_sample}.bam"
    params:
        model=DORADO_MODEL
    log:
        "logs/dorado_basecall/{raw_sample}.log"
    benchmark:
        "benchmarks/phase0/dorado_basecall/{raw_sample}.tsv"
    threads: 4
    resources:
        gpu=1
    shell:
        """
        workflow/scripts/Dorado_basecall.sh \
            -i {input.pod5} \
            -o {output.bam} \
            -m {params.model} \
            -t {threads} \
            2>&1 | tee {log}
        """


rule uncalled4_align:
    """Align raw nanopore signals to pore model using Uncalled4 (BAM output, all strands).
    Used by perbase_error (strand filtering happens inside that script via samtools view).
    --bam-in requires a reference-aligned BAM that retains the dorado mv/ts/pi/sp/ns tags.
    Map_reads.sh preserves these tags via samtools fastq -T + minimap2 -y, so the tags
    survive through aligned_reads → filtered_alignments.
    Command syntax: uncalled4 align --bam-in --ref --reads -o (confirmed from js4004).
    """
    input:
        bam="data/filtered_alignments/{genome}/{raw_sample}.bam",
        pod5=lambda wildcards: pod5_input(wildcards.raw_sample),
        genome="resources/genomes/{genome}.fa"
    output:
        bam="data/uncalled4/{genome}/{raw_sample}.bam",
        bai="data/uncalled4/{genome}/{raw_sample}.bam.bai"
    log:
        "logs/uncalled4_align/{genome}/{raw_sample}.log"
    benchmark:
        "benchmarks/phase0/uncalled4_align/{genome}/{raw_sample}.tsv"
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


rule uncalled4_convert_tsv:
    """Convert Uncalled4 BAM to a strand-specific DTW TSV using uncalled4 convert.
    Avoids re-running signal alignment — all DTW data is already in the BAM tags.
    Pre-filters to one strand with samtools view before converting.
    The reference is passed explicitly as --ref (via -g to the wrapper script).
    uncalled4 needs it for the pore model k-mer behind dtw.model_diff, and its
    fallback is the path recorded in the BAM header at alignment time, which does
    not resolve for a BAM produced by a different run or on a different machine.
    Note the flag is --ref, not the --ref-index named in uncalled4's own error.
    dtw.model_diff = observed - model current, in normalized (not pA) units;
    consumed by perbase_signal_deviation.py.
    The BAM comes from uncalled4_align for a pod5 sample, or straight from the
    ^u path for a sample supplied already signal-aligned.
    """
    input:
        bam=uncalled4_bam,
        bai=uncalled4_bai,
        genome="resources/genomes/{genome}.fa"
    output:
        tsv="data/uncalled4_tsv/{genome}/{raw_sample}_{strand}.tsv"
    log:
        "logs/uncalled4_convert_tsv/{genome}/{raw_sample}_{strand}.log"
    benchmark:
        "benchmarks/phase0/uncalled4_convert_tsv/{genome}/{raw_sample}_{strand}.tsv"
    wildcard_constraints:
        strand="for|rev"
    threads: 4
    shell:
        """
        workflow/scripts/Uncalled4_convert_tsv.sh \
            -i {input.bam} \
            -o {output.tsv} \
            -s {wildcards.strand} \
            -g {input.genome} \
            -t {threads} \
            2>&1 | tee {log}
        """
