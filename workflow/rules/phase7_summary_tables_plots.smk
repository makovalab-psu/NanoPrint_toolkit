# Phase 7: Summary Tables & Plots
# Rules for generating summary tables and visualizations


def get_correlation_files_for_summary(wildcards):
    """Get all pairwise correlation files (both strands) for a given genome."""
    return expand(
        "tables/correlation/{genome}/{raw_sample_a}_vs_{raw_sample_b}_{strand}.txt",
        genome=wildcards.genome,
        raw_sample_a=RAW_SAMPLES,
        raw_sample_b=RAW_SAMPLES,
        strand=["for", "rev"]
    )


rule correlation:
    """Correlate per-base error between two samples via random subsampling."""
    input:
        file_a="data/perbase_error/{genome}/{raw_sample_a}_{strand}.txt.gz",
        file_b="data/perbase_error/{genome}/{raw_sample_b}_{strand}.txt.gz",
        fai="resources/genomes/{genome}.fa.fai"
    output:
        correlation="tables/correlation/{genome}/{raw_sample_a}_vs_{raw_sample_b}_{strand}.txt"
    params:
        subsamples=config.get("correlation_subsamples", 10000)
    log:
        "logs/correlation/{genome}/{raw_sample_a}_vs_{raw_sample_b}_{strand}.log"
    benchmark:
        "benchmarks/phase7/correlation/{genome}/{raw_sample_a}_vs_{raw_sample_b}_{strand}.tsv"
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


rule summarize_correlation:
    """Compute pairwise Spearman/Pearson correlations and generate heatmap plots."""
    input:
        get_correlation_files_for_summary
    output:
        table="tables/perbase_error_correlation/{genome}/Pairwise_correlation_table.csv",
        pdf="plots/perbase_error_correlation/{genome}/Pairwise_correlation_heatmap.pdf"
    log:
        "logs/summarize_correlation/{genome}.log"
    benchmark:
        "benchmarks/phase7/summarize_correlation/{genome}.tsv"
    shell:
        """
        Rscript workflow/scripts/Summarize_correlation.R \
            {output.table} \
            {output.pdf} \
            {input} \
            2>&1 | tee {log}
        """


rule summary_histogram_plot:
    """Generate summary histogram PDF (Read Length, MAPQ, Insertion, Deletion, Coverage)."""
    input:
        txt="data/{alignment}/{genome}/{raw_sample}_histograms.txt"
    output:
        pdf="plots/histograms/{alignment}/{genome}/{raw_sample}_histograms.pdf"
    log:
        "logs/histograms_plot/{alignment}/{genome}/{raw_sample}.log"
    benchmark:
        "benchmarks/phase7/summary_histogram_plot/{alignment}/{genome}/{raw_sample}.tsv"
    wildcard_constraints:
        alignment="aligned_reads|filtered_alignments"
    shell:
        """
        Rscript workflow/scripts/Plot_summary_histograms.R \
            {input.txt} \
            {output.pdf} \
            2>&1 | tee {log}
        """


rule plot_annotation:
    """Plot mean coverage, per-base error, and reactivity around genomic features."""
    input:
        for_strand="data/annotations_averaged/{genome}/{feature}/{sample}_for.txt.gz",
        rev_strand="data/annotations_averaged/{genome}/{feature}/{sample}_rev.txt.gz"
    output:
        pdf="plots/annotations_averaged/{genome}/{feature}/{sample}.pdf"
    log:
        "logs/plot_annotation/{genome}/{feature}/{sample}.log"
    benchmark:
        "benchmarks/phase7/plot_annotation/{genome}/{feature}/{sample}.tsv"
    shell:
        """
        Rscript workflow/scripts/Plot_annotation.R \
            {input.for_strand} \
            {input.rev_strand} \
            {output.pdf} \
            2>&1 | tee {log}
        """


rule alignment_stats_table:
    """Combine per-sample alignment statistics into a single CSV table (two rows per sample: raw and filtered)."""
    input:
        expand("tables/alignment_stats/{genome}/{raw_sample}.txt",
               genome=GENOMES, raw_sample=RAW_SAMPLES)
    output:
        csv="tables/alignment_stats_table.csv"
    log:
        "logs/alignment_stats_table.log"
    benchmark:
        "benchmarks/phase7/alignment_stats_table.tsv"
    shell:
        """
        workflow/scripts/alignment_stats_table.sh \
            -o {output.csv} \
            {input} \
            2>&1 | tee {log}
        """


rule read_stats_table:
    """Combine per-sample read statistics into a single CSV table."""
    input:
        expand("tables/read_stats/{raw_sample}.txt", raw_sample=RAW_SAMPLES)
    output:
        csv="tables/read_stats_table.csv"
    log:
        "logs/read_stats_table.log"
    benchmark:
        "benchmarks/phase7/read_stats_table.tsv"
    shell:
        """
        workflow/scripts/read_stats_table.sh \
            -o {output.csv} \
            {input} \
            2>&1 | tee {log}
        """
