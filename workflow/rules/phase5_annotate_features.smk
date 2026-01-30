# Phase 5: Feature Annotation
# Rules for calculating signal around genomic features


rule split_features_by_chr:
    """Split feature BED file by chromosome for parallelization."""
    input:
        bed="resources/features/{feature}.bed"
    output:
        done="resources/features/{feature}_by_chr/.done"
    params:
        outdir="resources/features/{feature}_by_chr"
    log:
        "logs/split_features/{feature}.log"
    shell:
        """
        mkdir -p {params.outdir}
        workflow/scripts/Split_by_chr.sh \
            -i {input.bed} \
            2>&1 | tee {log}
        # Move split files to output directory
        mv resources/features/{wildcards.feature}_*.bed {params.outdir}/
        touch {output.done}
        """


rule feature_chr_file:
    """Declare individual chromosome BED files produced by split_features_by_chr."""
    input:
        done="resources/features/{feature}_by_chr/.done"
    output:
        file="resources/features/{feature}_by_chr/{feature}_{chr}.bed"
    shell:
        """
        # File was created by split_features_by_chr, just verify it exists
        test -f {output.file}
        """


def get_feature_chromosomes(wildcards):
    """Get list of chromosomes for a feature."""
    return FEATURE_CHROMOSOMES[wildcards.feature]


def get_annotation_chr_files(wildcards):
    """Get all chromosome annotation files for merging."""
    chrs = FEATURE_CHROMOSOMES[wildcards.feature]
    return expand(
        "data/annotations/{genome}/{feature}/{sample}_{strand}_{chr}.txt.gz",
        genome=wildcards.genome,
        feature=wildcards.feature,
        sample=wildcards.sample,
        strand=wildcards.strand,
        chr=chrs
    )


def get_annotate_inputs(wildcards):
    """Get input files for annotate_features rule."""
    treatment = get_treatment(wildcards.sample)
    control = get_control(wildcards.sample)
    return {
        "bed": f"resources/features/{wildcards.feature}_by_chr/{wildcards.feature}_{wildcards.chr}.bed",
        "treatment": f"data/perbase_error_by_chr/{wildcards.genome}/{treatment}_{wildcards.strand}/{treatment}_{wildcards.strand}_{wildcards.chr}.txt",
        "control": f"data/perbase_error_by_chr/{wildcards.genome}/{control}_{wildcards.strand}/{control}_{wildcards.strand}_{wildcards.chr}.txt",
        "reactivity": f"data/reactivity/{wildcards.genome}/{wildcards.sample}_{wildcards.strand}_{wildcards.chr}.txt.gz"
    }


rule annotate_features:
    """Calculate signal around genomic features (per chromosome)."""
    input:
        unpack(get_annotate_inputs)
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
