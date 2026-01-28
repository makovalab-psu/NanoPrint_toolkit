# Phase 4: Output Formats & Analysis
# Rules for generating bedGraph, bigWig, and density analysis


rule reactivity_to_bedgraph:
    """Convert reactivity data to bedGraph format."""
    input:
        reactivity="data/reactivity/{genome}/{sample}_{strand}_{chr}.txt.gz"
    output:
        bedgraph="data/bg/{genome}/{sample}_{strand}_{chr}.bg"
    log:
        "logs/bedgraph/{genome}/{sample}_{strand}_{chr}.log"
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
        bedgraph="data/bg_by_chr/{genome}/{sample}_{strand}/{sample}_{strand}_{chr}.bg",
        fai="resources/genomes/{genome}.fa.fai"
    output:
        density="data/windows/{genome}/window_size_{size}/significance_threshold_{sig}/{sample}_{strand}_{chr}.bg"
    log:
        "logs/density/{genome}/{sample}_{strand}_{chr}_w{size}_s{sig}.log"
    wildcard_constraints:
        strand="for|rev"
    shell:
        """
        workflow/scripts/react_dens.sh \
            -i {input.bedgraph} \
            -f {input.fai} \
            -w {wildcards.size} \
            -s {wildcards.sig} \
            -o {output.density} \
            2>&1 | tee {log}
        """

rule merge_density:
    """Merge chromosome-split density files in genome order."""
    input:
        files=expand("data/windows/{genome}/window_size_{size}/significance_threshold_{sig}/{sample}_{strand}_{chr}.bg", chr = CHR)
        fai="resources/genomes/{genome}.fa.fai"
    output:
        merged="data/windows_merged/{genome}/window_size_{size}/significance_threshold_{sig}/{sample}_{strand}.bg"
    log:
        "logs/merge_density/{genome}/{sample}_{strand}_w{size}_s{sig}.log"
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
        bedgraph="data/bg_by_chr/{genome}/{sample}_{strand}/{sample}_{strand}_{chr}.bg",
        fai="resources/genomes/{genome}.fa.fai"
    output:
        bigwig="data/bw/{genome}/significance_threshold_{sig}/{sample}_{strand}_{chr}.bw"
    log:
        "logs/bigwig/{genome}/{sample}_{strand}_{chr}_s{sig}.log"
    wildcard_constraints:
        strand="for|rev"
    shell:
        """
        workflow/scripts/bg_to_bw.sh \
            -i {input.bedgraph} \
            -f {input.fai} \
            -s {wildcards.sig} \
            -o {output.bigwig} \
            2>&1 | tee {log}
        """

rule merge_bigwig:
    """Merge chromosome-split bigWig files in genome order."""
    input:
        files=expand("data/bw/{genome}/significance_threshold_{sig}/{sample}_{strand}_{chr}.bw", chr = CHR)
        fai="resources/genomes/{genome}.fa.fai"
    output:
        merged="data/bw_merged/{genome}/significance_threshold_{sig}/{sample}_{strand}.bw"
    log:
        "logs/merge_bigwig/{genome}/{sample}_{strand}_s{sig}.log"
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
