# Phase 3: Reactivity Calculation
# Rules for calculating reactivity (treatment - control) per chromosome


def get_reactivity_inputs(wildcards):
    """Get input files for calculate_reactivity rule."""
    treatment = get_treatment(wildcards.sample)
    control = get_control(wildcards.sample)
    return {
        "treatment": f"data/perbase_error_by_chr/{wildcards.genome}/{treatment}_{wildcards.strand}/{treatment}_{wildcards.strand}_{wildcards.chr}.txt.gz",
        "control": f"data/perbase_error_by_chr/{wildcards.genome}/{control}_{wildcards.strand}/{control}_{wildcards.strand}_{wildcards.chr}.txt.gz"
    }


rule calculate_reactivity:
    """Calculate reactivity by comparing treatment to control samples (per chromosome)."""
    input:
        unpack(get_reactivity_inputs)
    output:
        reactivity=wrap_output("reactivity", "data/reactivity/{genome}/{sample}_{strand}_{chr}.txt.gz")
    log:
        "logs/reactivity/{genome}/{sample}_{strand}_{chr}.log"
    benchmark:
        "benchmarks/phase3/calculate_reactivity/{genome}/{sample}_{strand}_{chr}.tsv"
    wildcard_constraints:
        strand="for|rev"
    shell:
        """
        workflow/scripts/Calculate_reactivity.sh \
            -p {input.treatment} \
            -m {input.control} \
            -o {output.reactivity} \
            2>&1 | tee {log}
        """
