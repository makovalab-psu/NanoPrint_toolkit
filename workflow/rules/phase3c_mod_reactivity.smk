# Phase 3c: Modification-Rate Reactivity
# Treatment minus control modified fraction, both scored against the same k-mer
# model (phase 2c): the control's own ({model_set} = matched) or a ^mod-pool's.
# The control term is the false-positive rate of the test on unmodified DNA,
# measured on reads the model never saw.
#
# Whole-genome files, merged by contig and position in mod_reactivity.py, so there
# is no per-chromosome split and no genome-specific rule.


def get_mod_reactivity_inputs(wildcards):
    """Treatment and control calls for a relationship, both against one model."""
    treatment = get_treatment(wildcards.sample)
    control = get_control(wildcards.sample)
    model = control if wildcards.model_set == "matched" else wildcards.model_set
    base = f"data/perbase_mod/{wildcards.backend}/{wildcards.genome}/{wildcards.thr}"
    return {
        "treatment": f"{base}/{treatment}_vs_{model}_{wildcards.strand}.txt.gz",
        "control": f"{base}/{control}_vs_{model}_{wildcards.strand}.txt.gz",
        "fai": f"resources/genomes/{wildcards.genome}.fa.fai",
    }


rule calculate_mod_reactivity:
    """Calculate modification-rate reactivity: treatment - control modified fraction."""
    input:
        unpack(get_mod_reactivity_inputs)
    output:
        reactivity="data/mod_reactivity/{backend}/{genome}/{thr}/{model_set}/{sample}_{strand}.txt.gz"
    params:
        # Same key and meaning as phases 3 and 3b: per-strand coverage, both samples.
        cov=config.get("reactivity_cov_threshold", 10)
    log:
        "logs/mod_reactivity/{backend}/{genome}/{thr}/{model_set}/{sample}_{strand}.log"
    benchmark:
        "benchmarks/phase3c/calculate_mod_reactivity/{backend}/{genome}/{thr}/{model_set}/{sample}_{strand}.tsv"
    wildcard_constraints:
        backend="uncalled4|remora",
        strand="for|rev",
        thr=MOD_THR_CONSTRAINT,
        model_set="|".join(["matched"] + [re.escape(m) for m in MOD_POOLS])
    shell:
        """
        python3 workflow/scripts/mod_reactivity.py \
            -p {input.treatment} \
            -m {input.control} \
            -g {input.fai} \
            -c {params.cov} \
            -o {output.reactivity} \
            2>&1 | tee {log}
        """
