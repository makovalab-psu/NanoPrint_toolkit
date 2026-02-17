# Phase 4: Output Formats & Analysis
# Rules for generating bedGraph, bigWig, and density analysis


def get_density_chr_files(wildcards):
    """Get all chromosome density files for merging."""
    chrs = CHROMOSOMES[wildcards.genome]
    return expand(
        "data/windows/{genome}/window_size_{size}/significance_threshold_{sig}/{sample}_{strand}_{chr}.bg",
        genome=wildcards.genome,
        size=wildcards.size,
        sig=wildcards.sig,
        sample=wildcards.sample,
        strand=wildcards.strand,
        chr=chrs
    )


def get_bigwig_chr_files(wildcards):
    """Get all chromosome bigWig files for merging."""
    chrs = CHROMOSOMES[wildcards.genome]
    return expand(
        "data/bw/{genome}/significance_threshold_{sig}/{sample}_{strand}_{chr}.bw",
        genome=wildcards.genome,
        sig=wildcards.sig,
        sample=wildcards.sample,
        strand=wildcards.strand,
        chr=chrs
    )


rule reactivity_to_bedgraph:
    """Convert reactivity data to bedGraph format."""
    input:
        reactivity="data/reactivity/{genome}/{sample}_{strand}_{chr}.txt.gz"
    output:
        bedgraph=wrap_output("bg", "data/bg/{genome}/{sample}_{strand}_{chr}.bg")
    log:
        "logs/bedgraph/{genome}/{sample}_{strand}_{chr}.log"
    benchmark:
        "benchmarks/phase4/reactivity_to_bedgraph/{genome}/{sample}_{strand}_{chr}.tsv"
    wildcard_constraints:
        strand="for|rev"
    shell:
        """
        workflow/scripts/react_to_bg.sh \
            -i {input.reactivity} \
            -o {output.bedgraph} \
            2>&1 | tee {log}
        """

rule reactivity_density:
    """Calculate reactive nucleotide density in genomic windows (per chromosome)."""
    input:
        bedgraph="data/bg/{genome}/{sample}_{strand}_{chr}.bg",
        fai="resources/genomes/{genome}.fa.fai"
    output:
        density=wrap_output("windows", "data/windows/{genome}/window_size_{size}/significance_threshold_{sig}/{sample}_{strand}_{chr}.bg")
    log:
        "logs/density/{genome}/{sample}_{strand}_{chr}_w{size}_s{sig}.log"
    benchmark:
        "benchmarks/phase4/reactivity_density/{genome}/{sample}_{strand}_{chr}_{size}_{sig}.tsv"
    wildcard_constraints:
        strand="for|rev"
    shell:
        """
        workflow/scripts/react_dens.sh \
            -i {input.bedgraph} \
            -g {input.fai} \
            -w {wildcards.size} \
            -p {wildcards.sig} \
            -o {output.density} \
            2>&1 | tee {log}
        """

rule merge_density:
    """Merge chromosome-split density files in genome order."""
    input:
        files=get_density_chr_files,
        fai="resources/genomes/{genome}.fa.fai"
    output:
        merged="data/windows_merged/{genome}/window_size_{size}/significance_threshold_{sig}/{sample}_{strand}.bg"
    log:
        "logs/merge_density/{genome}/{sample}_{strand}_w{size}_s{sig}.log"
    benchmark:
        "benchmarks/phase4/merge_density/{genome}/{sample}_{strand}_{size}_{sig}.tsv"
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

rule bedgraph_to_bigwig:
    """Convert bedGraph to bigWig format (per chromosome)."""
    input:
        bedgraph="data/bg/{genome}/{sample}_{strand}_{chr}.bg",
        fai="resources/genomes/{genome}.fa.fai"
    output:
        bigwig=wrap_output("bw", "data/bw/{genome}/significance_threshold_{sig}/{sample}_{strand}_{chr}.bw")
    log:
        "logs/bigwig/{genome}/{sample}_{strand}_{chr}_s{sig}.log"
    benchmark:
        "benchmarks/phase4/bedgraph_to_bigwig/{genome}/{sample}_{strand}_{chr}_{sig}.tsv"
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

rule merge_bigwig:
    """Merge chromosome-split bigWig files in genome order."""
    input:
        files=get_bigwig_chr_files,
        fai="resources/genomes/{genome}.fa.fai"
    output:
        merged="data/bw_merged/{genome}/significance_threshold_{sig}/{sample}_{strand}.bw"
    log:
        "logs/merge_bigwig/{genome}/{sample}_{strand}_s{sig}.log"
    benchmark:
        "benchmarks/phase4/merge_bigwig/{genome}/{sample}_{strand}_{sig}.tsv"
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
# Mean reactivity in genomic windows
# ============================================================================

def get_mean_bg_chr_files(wildcards):
    """Get all chromosome mean bedGraph files for merging."""
    chrs = CHROMOSOMES[wildcards.genome]
    return expand(
        "data/bg_mean/{genome}/window_size_{mean_size}/{sample}_{strand}_{chr}.bg",
        genome=wildcards.genome,
        mean_size=wildcards.mean_size,
        sample=wildcards.sample,
        strand=wildcards.strand,
        chr=chrs
    )


rule mean_reactivity_bedgraph:
    """Calculate mean reactivity in genomic windows (per chromosome)."""
    input:
        reactivity="data/reactivity/{genome}/{sample}_{strand}_{chr}.txt.gz",
        fai="resources/genomes/{genome}.fa.fai"
    output:
        bedgraph=wrap_output("bg_mean", "data/bg_mean/{genome}/window_size_{mean_size}/{sample}_{strand}_{chr}.bg")
    log:
        "logs/mean_bedgraph/{genome}/{sample}_{strand}_{chr}_w{mean_size}.log"
    benchmark:
        "benchmarks/phase4/mean_reactivity_bedgraph/{genome}/{sample}_{strand}_{chr}_{mean_size}.tsv"
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

rule merge_mean_bedgraph:
    """Merge chromosome-split mean reactivity bedGraph files in genome order."""
    input:
        files=get_mean_bg_chr_files,
        fai="resources/genomes/{genome}.fa.fai"
    output:
        merged="data/bg_mean_merged/{genome}/window_size_{mean_size}/{sample}_{strand}.bg"
    log:
        "logs/merge_mean_bedgraph/{genome}/{sample}_{strand}_w{mean_size}.log"
    benchmark:
        "benchmarks/phase4/merge_mean_bedgraph/{genome}/{sample}_{strand}_{mean_size}.tsv"
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

rule mean_bedgraph_to_bigwig:
    """Convert merged mean reactivity bedGraph to bigWig format."""
    input:
        bedgraph="data/bg_mean_merged/{genome}/window_size_{mean_size}/{sample}_{strand}.bg",
        fai="resources/genomes/{genome}.fa.fai"
    output:
        bigwig="data/bw_mean_merged/{genome}/window_size_{mean_size}/{sample}_{strand}.bw"
    log:
        "logs/mean_bigwig/{genome}/{sample}_{strand}_w{mean_size}.log"
    benchmark:
        "benchmarks/phase4/mean_bedgraph_to_bigwig/{genome}/{sample}_{strand}_{mean_size}.tsv"
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
