# Phase 1: Read Mapping & Quality Control
# Rules for mapping, filtering, and QC statistics

import os
import glob as pyglob


def find_raw_reads(wildcards):
    """Find raw read file regardless of extension (.fastq, .fastq.gz, .bam)."""
    raw_sample = wildcards.raw_sample
    base_path = f"raw_data/{raw_sample}"

    # Check for each supported extension in order of preference
    for ext in [".fastq.gz", ".fastq", ".bam"]:
        if os.path.exists(base_path + ext):
            return base_path + ext

    # If no file found, return expected path (will fail with clear error)
    raise FileNotFoundError(
        f"No raw reads found for {raw_sample}. "
        f"Expected one of: {base_path}.fastq.gz, {base_path}.fastq, {base_path}.bam"
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
    """Extract histogram data and generate PDF plots from stats files."""
    input:
        stats="data/{alignment}/{genome}/{raw_sample}_stats.txt"
    output:
        txt="data/{alignment}/{genome}/{raw_sample}_histograms.txt",
        pdf="data/{alignment}/{genome}/{raw_sample}_histograms.pdf"
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
            -p {output.pdf} \
            2>&1 | tee {log}
        """
