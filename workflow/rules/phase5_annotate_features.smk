# Phase 5: Feature Annotation
# Rules for calculating signal around genomic features


checkpoint split_features_by_chr:
    """Split feature BED file by chromosome for parallelization."""
    input:
        bed="resources/features/{feature}.bed"
    output:
        directory("resources/features/{feature}_by_chr")
    log:
        "logs/split_features/{feature}.log"
    shell:
        """
        mkdir -p {output}
        workflow/scripts/Split_by_chr.sh \
            -i {input.bed} \
            2>&1 | tee {log}
        # Move split files to output directory
        mv resources/features/{wildcards.feature}_*.bed {output}/
        """


def get_feature_chromosomes(wildcards):
    """Get list of chromosomes from split feature checkpoint."""
    checkpoint_output = checkpoints.split_features_by_chr.get(
        feature=wildcards.feature
    ).output[0]
    import glob
    import os
    pattern = os.path.join(checkpoint_output, f"{wildcards.feature}_*.bed")
    files = glob.glob(pattern)
    chrs = [os.path.basename(f).replace(f"{wildcards.feature}_", "").replace(".bed", "") for f in files]
    return chrs


def get_annotation_chr_files(wildcards):
    """Get all chromosome annotation files for merging."""
    chrs = get_feature_chromosomes(wildcards)
    return expand(
        "data/annotations/{genome}/{feature}/{sample}_{strand}_{chr}.txt.gz",
        genome=wildcards.genome,
        feature=wildcards.feature,
        sample=wildcards.sample,
        strand=wildcards.strand,
        chr=chrs
    )


rule annotate_features:
    """Calculate signal around genomic features (per chromosome)."""
    input:
        bed="resources/features/{feature}_by_chr/{feature}_{chr}.bed",
        treatment="data/perbase_error_by_chr/{genome}/{treatment_sample}_{strand}/{treatment_sample}_{strand}_{chr}.txt",
        control="data/perbase_error_by_chr/{genome}/{control_sample}_{strand}/{control_sample}_{strand}_{chr}.txt",
        reactivity="data/reactivity/{genome}/{sample}_{strand}_{chr}.txt.gz"
    output:
        annotation="data/annotations/{genome}/{feature}/{sample}_{strand}_{chr}.txt.gz"
    params:
        n_windows=lambda wildcards: FEATURE_PARAMS.get(wildcards.feature, (1000, 10))[0],
        window_size=lambda wildcards: FEATURE_PARAMS.get(wildcards.feature, (1000, 10))[1]
    log:
        "logs/annotate_features/{genome}/{feature}/{sample}_{strand}_{chr}.log"
    wildcard_constraints:
        strand="for|rev"
    shell:
        """
        workflow/scripts/annotate_features.sh \
            -b {input.bed} \
            -t {input.treatment} \
            -c {input.control} \
            -r {input.reactivity} \
            -o {output.annotation} \
            -n {params.n_windows} \
            -w {params.window_size} \
            2>&1 | tee {log}
        """


rule merge_annotations:
    """Merge chromosome-split annotation files."""
    input:
        files=get_annotation_chr_files
    output:
        merged="data/annotations_merged/{genome}/{feature}/{sample}_{strand}.txt.gz"
    log:
        "logs/merge_annotations/{genome}/{feature}/{sample}_{strand}.log"
    wildcard_constraints:
        strand="for|rev"
    shell:
        """
        workflow/scripts/Merge_annotations.sh \
            -o {output.merged} \
            {input.files} \
            2>&1 | tee {log}
        """


rule average_annotations:
    """Average annotations by distance, sample, and strand."""
    input:
        merged="data/annotations_merged/{genome}/{feature}/{sample}_{strand}.txt.gz"
    output:
        averaged="data/annotations_averaged/{genome}/{feature}/{sample}_{strand}.txt.gz"
    log:
        "logs/average_annotations/{genome}/{feature}/{sample}_{strand}.log"
    wildcard_constraints:
        strand="for|rev"
    shell:
        """
        workflow/scripts/average_feature_annotation.sh \
            -i {input.merged} \
            -o {output.averaged} \
            2>&1 | tee {log}
        """
