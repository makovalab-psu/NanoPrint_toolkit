# Phase 2: Per-Base Error Quantification
# Rules for calculating per-base error rates by strand and correlation


def find_perbase_bam(wildcards):
    """Use the Uncalled4 BAM when there is one; the filtered BAM otherwise.

    Three cases, in the order they are tested:
      ^u sample   the supplied Uncalled4 BAM, read where it lies
      pod5 sample the Uncalled4 BAM that uncalled4_align produces
      otherwise   data/filtered_alignments, with no signal data in play
    """
    if has_signal(wildcards.raw_sample):
        return uncalled4_bam(wildcards)
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
