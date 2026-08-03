# Phase 2: Per-Base Error Quantification
# Rules for calculating per-base error rates by strand and correlation


def find_perbase_bam(wildcards):
    """Use Uncalled4 BAM when pod5 input is available; filtered BAM otherwise."""
    if has_pod5(wildcards.raw_sample):
        return f"data/uncalled4/{wildcards.genome}/{wildcards.raw_sample}.bam"
    return f"data/filtered_alignments/{wildcards.genome}/{wildcards.raw_sample}.bam"


rule perbase_error:
    """Calculate per-base error rates from filtered alignments or Uncalled4 BAM."""
    input:
        bam=find_perbase_bam,
        genome="resources/genomes/{genome}.fa"
    output:
        error="data/perbase_error/{genome}/{raw_sample}_{strand}.txt.gz"
    params:
        strand_flag=lambda wildcards: "--for" if wildcards.strand == "for" else "--rev"
    log:
        "logs/perbase_error/{genome}/{raw_sample}_{strand}.log"
    benchmark:
        "benchmarks/phase2/perbase_error/{genome}/{raw_sample}_{strand}.tsv"
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


# Note: split_perbase_by_chr rules are generated per-genome in genome_specific_rules.smk
# Note: correlation rule moved to phase7_summary_tables_plots.smk
