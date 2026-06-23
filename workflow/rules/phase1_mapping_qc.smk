# Phase 1: Read Mapping & Quality Control
# Rules for mapping, filtering, and QC statistics

import os
import glob as pyglob


def find_raw_reads(wildcards):
    """Return the raw read file for a sample, using the path from RAW_PATHS.

    For pod5 input: dorado_basecall (phase 0) produces data/basecalled/{sample}.bam,
    which is returned here so map_reads depends on it and triggers basecalling.
    For all other input types (fastq.gz, fastq, bam): the path is used directly.
    """
    raw_sample = wildcards.raw_sample

    # Pod5 input: basecalling handled by dorado_basecall rule (phase 0)
    if has_pod5(raw_sample):
        return f"data/basecalled/{raw_sample}.bam"

    path = RAW_PATHS.get(raw_sample)
    if path and os.path.exists(path):
        return path

    raise FileNotFoundError(
        f"No raw reads found for '{raw_sample}'. "
        f"Configured path: {path!r}. "
        f"Check the ^r line in CONFIG and ensure the file exists."
    )


rule genome_faidx:
    """Create FASTA index file for reference genome."""
    input:
        fa="resources/genomes/{genome}.fa"
    output:
        fai="resources/genomes/{genome}.fa.fai"
    log:
        "logs/genome_faidx/{genome}.log"
    benchmark:
        "benchmarks/phase1/genome_faidx/{genome}.tsv"
    shell:
        """
        samtools faidx {input.fa} 2>&1 | tee {log}
        """


rule read_stats:
    """Calculate sequencing statistics from raw reads."""
    input:
        reads=find_raw_reads
    output:
        stats="tables/read_stats/{raw_sample}.txt"
    log:
        "logs/read_stats/{raw_sample}.log"
    benchmark:
        "benchmarks/phase1/read_stats/{raw_sample}.tsv"
    shell:
        """
        workflow/scripts/Read_stats.sh \
            -i {input.reads} \
            -o {output.stats} \
            2>&1 | tee {log}
        """


rule map_reads:
    """Map raw reads to reference genome using minimap2."""
    input:
        reads=find_raw_reads,
        genome="resources/genomes/{genome}.fa"
    output:
        bam=wrap_output("aligned_reads_bam", "data/aligned_reads/{genome}/{raw_sample}.bam")
    log:
        "logs/map_reads/{genome}/{raw_sample}.log"
    benchmark:
        "benchmarks/phase1/map_reads/{genome}/{raw_sample}.tsv"
    threads: 8
    shell:
        """
        workflow/scripts/Map_reads.sh \
            -i {input.reads} \
            -g {input.genome} \
            -o {output.bam} \
            -t {threads} \
            2>&1 | tee {log}
        """


rule filter_alignments:
    """Filter alignments by mapping quality and remove secondary/supplementary."""
    input:
        bam="data/aligned_reads/{genome}/{raw_sample}.bam"
    output:
        bam="data/filtered_alignments/{genome}/{raw_sample}.bam"
    log:
        "logs/filter_alignments/{genome}/{raw_sample}.log"
    benchmark:
        "benchmarks/phase1/filter_alignments/{genome}/{raw_sample}.tsv"
    shell:
        """
        workflow/scripts/Filter_alignments.sh \
            -i {input.bam} \
            -o {output.bam} \
            2>&1 | tee {log}
        """


rule alignment_stats:
    """Generate alignment statistics from raw and filtered BAMs."""
    input:
        raw_bam="data/aligned_reads/{genome}/{raw_sample}.bam",
        filtered_bam="data/filtered_alignments/{genome}/{raw_sample}.bam"
    output:
        raw_flagstats="data/aligned_reads/{genome}/{raw_sample}_flagstats.txt",
        raw_stats="data/aligned_reads/{genome}/{raw_sample}_stats.txt",
        filtered_flagstats="data/filtered_alignments/{genome}/{raw_sample}_flagstats.txt",
        filtered_stats="data/filtered_alignments/{genome}/{raw_sample}_stats.txt",
        table="tables/alignment_stats/{genome}/{raw_sample}.txt"
    log:
        "logs/alignment_stats/{genome}/{raw_sample}.log"
    benchmark:
        "benchmarks/phase1/alignment_stats/{genome}/{raw_sample}.tsv"
    shell:
        """
        workflow/scripts/Alignment_stats.sh \
            -a {input.raw_bam} \
            -f {input.filtered_bam} \
            -o {output.table} \
            2>&1 | tee {log}
        """


rule histograms:
    """Extract histogram data from samtools stats files."""
    input:
        stats="data/{alignment}/{genome}/{raw_sample}_stats.txt"
    output:
        txt="data/{alignment}/{genome}/{raw_sample}_histograms.txt"
    log:
        "logs/histograms/{alignment}/{genome}/{raw_sample}.log"
    benchmark:
        "benchmarks/phase1/histograms/{alignment}/{genome}/{raw_sample}.tsv"
    wildcard_constraints:
        alignment="aligned_reads|filtered_alignments"
    shell:
        """
        workflow/scripts/Make_histograms.sh \
            -i {input.stats} \
            -o {output.txt} \
            2>&1 | tee {log}
        """
