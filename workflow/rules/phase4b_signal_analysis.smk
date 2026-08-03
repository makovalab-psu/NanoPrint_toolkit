# Phase 4b: Signal Reactivity Output Formats
# bedGraph and bigWig generation for signal reactivity data.
# Reuses all existing phase 4 scripts (react_to_bg.sh, bg_to_bw.sh,
# Merge_bigwig.sh) with different input/output paths.


def get_signal_bigwig_chr_files(wildcards):
    """Get per-chromosome signal bigWig files for merging."""
    chrs = CHROMOSOMES[wildcards.genome]
    return expand(
        "data/signal_bw/{genome}/significance_threshold_{sig}/{sample}_{strand}_{chr}.bw",
        genome=wildcards.genome,
        sig=wildcards.sig,
        sample=wildcards.sample,
        strand=wildcards.strand,
        chr=chrs,
    )


def get_signal_mean_bg_chr_files(wildcards):
    """Get per-chromosome signal mean bedGraph files for merging."""
    chrs = CHROMOSOMES[wildcards.genome]
    return expand(
        "data/signal_bg_mean/{genome}/window_size_{mean_size}/{sample}_{strand}_{chr}.bg",
        genome=wildcards.genome,
        mean_size=wildcards.mean_size,
        sample=wildcards.sample,
        strand=wildcards.strand,
        chr=chrs,
    )


rule signal_reactivity_to_bedgraph:
    """Convert signal reactivity to bedGraph format."""
    input:
        reactivity="data/signal_reactivity/{genome}/{sample}_{strand}_{chr}.txt.gz"
    output:
        bedgraph=wrap_output(
            "signal_bg",
            "data/signal_bg/{genome}/{sample}_{strand}_{chr}.bg"
        )
    log:
        "logs/signal_bedgraph/{genome}/{sample}_{strand}_{chr}.log"
    benchmark:
        "benchmarks/phase4b/signal_reactivity_to_bedgraph/{genome}/{sample}_{strand}_{chr}.tsv"
    wildcard_constraints:
        strand="for|rev"
    shell:
        """
        workflow/scripts/react_to_bg.sh \
            -i {input.reactivity} \
            -o {output.bedgraph} \
            2>&1 | tee {log}
        """


rule signal_bedgraph_to_bigwig:
    """Convert signal bedGraph to bigWig format (per chromosome)."""
    input:
        bedgraph="data/signal_bg/{genome}/{sample}_{strand}_{chr}.bg",
        fai="resources/genomes/{genome}.fa.fai"
    output:
        bigwig=wrap_output(
            "signal_bw",
            "data/signal_bw/{genome}/significance_threshold_{sig}/{sample}_{strand}_{chr}.bw"
        )
    log:
        "logs/signal_bigwig/{genome}/{sample}_{strand}_{chr}_s{sig}.log"
    benchmark:
        "benchmarks/phase4b/signal_bedgraph_to_bigwig/{genome}/{sample}_{strand}_{chr}_{sig}.tsv"
    wildcard_constraints:
        strand="for|rev"
    shell:
        """
        workflow/scripts/bg_to_bw.sh \
            -i {input.bedgraph} \
            -g {input.fai} \
            -p {wildcards.sig} \
            -o {output.bigwig} \
            2>&1 | tee {log}
        """


rule merge_signal_bigwig:
    """Merge chromosome-split signal bigWig files in genome order."""
    input:
        files=get_signal_bigwig_chr_files,
        fai="resources/genomes/{genome}.fa.fai"
    output:
        merged="data/signal_bw_merged/{genome}/significance_threshold_{sig}/{sample}_{strand}.bw"
    log:
        "logs/merge_signal_bigwig/{genome}/{sample}_{strand}_s{sig}.log"
    benchmark:
        "benchmarks/phase4b/merge_signal_bigwig/{genome}/{sample}_{strand}_{sig}.tsv"
    wildcard_constraints:
        strand="for|rev"
    shell:
        """
        workflow/scripts/Merge_bigwig.sh \
            -g {input.fai} \
            -o {output.merged} \
            {input.files} \
            2>&1 | tee {log}
        """


# ============================================================================
# Mean signal deviation in genomic windows
# ============================================================================

rule mean_signal_reactivity_bedgraph:
    """Calculate mean signal reactivity in genomic windows (per chromosome)."""
    input:
        reactivity="data/signal_reactivity/{genome}/{sample}_{strand}_{chr}.txt.gz",
        fai="resources/genomes/{genome}.fa.fai"
    output:
        bedgraph=wrap_output(
            "signal_bg_mean",
            "data/signal_bg_mean/{genome}/window_size_{mean_size}/{sample}_{strand}_{chr}.bg"
        )
    log:
        "logs/mean_signal_bedgraph/{genome}/{sample}_{strand}_{chr}_w{mean_size}.log"
    benchmark:
        "benchmarks/phase4b/mean_signal_reactivity_bedgraph/{genome}/{sample}_{strand}_{chr}_{mean_size}.tsv"
    wildcard_constraints:
        strand="for|rev"
    shell:
        """
        workflow/scripts/react_mean_bg.sh \
            -i {input.reactivity} \
            -g {input.fai} \
            -w {wildcards.mean_size} \
            -o {output.bedgraph} \
            2>&1 | tee {log}
        """


rule merge_mean_signal_bedgraph:
    """Merge chromosome-split mean signal bedGraph files in genome order."""
    input:
        files=get_signal_mean_bg_chr_files,
        fai="resources/genomes/{genome}.fa.fai"
    output:
        merged="data/signal_bg_mean_merged/{genome}/window_size_{mean_size}/{sample}_{strand}.bg"
    log:
        "logs/merge_mean_signal_bedgraph/{genome}/{sample}_{strand}_w{mean_size}.log"
    benchmark:
        "benchmarks/phase4b/merge_mean_signal_bedgraph/{genome}/{sample}_{strand}_{mean_size}.tsv"
    wildcard_constraints:
        strand="for|rev"
    shell:
        """
        workflow/scripts/Merge_density.sh \
            -g {input.fai} \
            -o {output.merged} \
            {input.files} \
            2>&1 | tee {log}
        """


rule mean_signal_bedgraph_to_bigwig:
    """Convert merged mean signal bedGraph to bigWig format."""
    input:
        bedgraph="data/signal_bg_mean_merged/{genome}/window_size_{mean_size}/{sample}_{strand}.bg",
        fai="resources/genomes/{genome}.fa.fai"
    output:
        bigwig="data/signal_bw_mean_merged/{genome}/window_size_{mean_size}/{sample}_{strand}.bw"
    log:
        "logs/mean_signal_bigwig/{genome}/{sample}_{strand}_w{mean_size}.log"
    benchmark:
        "benchmarks/phase4b/mean_signal_bedgraph_to_bigwig/{genome}/{sample}_{strand}_{mean_size}.tsv"
    wildcard_constraints:
        strand="for|rev"
    shell:
        """
        workflow/scripts/mean_bg_to_bw.sh \
            -i {input.bedgraph} \
            -g {input.fai} \
            -o {output.bigwig} \
            2>&1 | tee {log}
        """
