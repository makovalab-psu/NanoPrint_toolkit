# Phase 2: Per-Base Error Quantification
# Rules for calculating per-base error rates by strand and correlation


rule perbase_error:
    """Calculate per-base error rates from filtered alignments."""
    input:
        bam="data/filtered_alignments/{genome}/{sample}.bam",
        genome="resources/genomes/{genome}.fa"
    output:
        error="data/perbase_error/{genome}/{sample}_{strand}.txt.gz"
    params:
        strand_flag=lambda wildcards: "--for" if wildcards.strand == "for" else "--rev"
    log:
        "logs/perbase_error/{genome}/{sample}_{strand}.log"
    wildcard_constraints:
        strand="for|rev"
    shell:
        """
        workflow/scripts/perbase_error.sh \
            {params.strand_flag} \
            -i {input.bam} \
            -g {input.genome} \
            -o {output.error} \
            2>&1 | tee {log}
        """


checkpoint split_perbase_by_chr:
    """Split per-base error file by chromosome for parallelization."""
    input:
        error="data/perbase_error/{genome}/{sample}_{strand}.txt.gz"
    output:
        directory("data/perbase_error_by_chr/{genome}/{sample}_{strand}")
    log:
        "logs/split_perbase_by_chr/{genome}/{sample}_{strand}.log"
    wildcard_constraints:
        strand="for|rev"
    shell:
        """
        mkdir -p {output}
        workflow/scripts/Split_by_chr.sh \
            -i {input.error} \
            2>&1 | tee {log}
        # Move split files to output directory
        mv data/perbase_error/{wildcards.genome}/{wildcards.sample}_{wildcards.strand}_*.txt {output}/
        """


rule correlation:
    """Correlate per-base error between two samples via random subsampling."""
    input:
        file_a="data/perbase_error/{genome}/{sample_a}_{strand}.txt.gz",
        file_b="data/perbase_error/{genome}/{sample_b}_{strand}.txt.gz",
        fai="resources/genomes/{genome}.fa.fai"
    output:
        correlation="tables/correlation/{genome}/{sample_a}_vs_{sample_b}_{strand}.txt"
    params:
        subsamples=config.get("correlation_subsamples", 10000)
    log:
        "logs/correlation/{genome}/{sample_a}_vs_{sample_b}_{strand}.log"
    wildcard_constraints:
        strand="for|rev"
    shell:
        """
        workflow/scripts/Correlation.sh \
            -a {input.file_a} \
            -b {input.file_b} \
            -g {input.fai} \
            -o {output.correlation} \
            -s {params.subsamples} \
            2>&1 | tee {log}
        """
