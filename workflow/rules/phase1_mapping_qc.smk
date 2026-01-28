# Phase 1: Read Mapping & Quality Control
# Rules for mapping, filtering, and QC statistics


rule read_stats:
    """Calculate sequencing statistics from raw reads."""
    input:
        reads="raw_data/{sample}.fastq.gz"
    output:
        stats="tables/read_stats/{sample}.txt"
    log:
        "logs/read_stats/{sample}.log"
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
        reads="raw_data/{sample}.fastq.gz",
        genome="resources/genomes/{genome}.fa"
    output:
        bam="data/aligned_reads/{genome}/{sample}.bam"
    log:
        "logs/map_reads/{genome}/{sample}.log"
    shell:
        """
        workflow/scripts/Map_reads.sh \
            -i {input.reads} \
            -g {input.genome} \
            -o {output.bam} \
            2>&1 | tee {log}
        """


rule filter_alignments:
    """Filter alignments by mapping quality and remove secondary/supplementary."""
    input:
        bam="data/aligned_reads/{genome}/{sample}.bam"
    output:
        bam="data/filtered_alignments/{genome}/{sample}.bam"
    log:
        "logs/filter_alignments/{genome}/{sample}.log"
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
        raw_bam="data/aligned_reads/{genome}/{sample}.bam",
        filtered_bam="data/filtered_alignments/{genome}/{sample}.bam"
    output:
        raw_flagstats="data/aligned_reads/{genome}/{sample}_flagstats.txt",
        raw_stats="data/aligned_reads/{genome}/{sample}_stats.txt",
        filtered_flagstats="data/filtered_alignments/{genome}/{sample}_flagstats.txt",
        filtered_stats="data/filtered_alignments/{genome}/{sample}_stats.txt",
        table="tables/alignment_stats/{genome}/{sample}.txt"
    log:
        "logs/alignment_stats/{genome}/{sample}.log"
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
        stats="data/{alignment}/{genome}/{sample}_stats.txt"
    output:
        txt="data/{alignment}/{genome}/{sample}_histograms.txt",
        pdf="data/{alignment}/{genome}/{sample}_histograms.pdf"
    log:
        "logs/histograms/{alignment}/{genome}/{sample}.log"
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
