# Phase 3: Reactivity Calculation
# Rules for calculating reactivity (treatment - control) per chromosome


def get_reactivity_inputs(wildcards):
    """Get input files for calculate_reactivity rule."""
    treatment = get_treatment(wildcards.sample)
    control = get_control(wildcards.sample)
    return {
        "treatment": f"data/perbase_error_by_chr/{wildcards.genome}/{treatment}_{wildcards.strand}/{treatment}_{wildcards.strand}_{wildcards.chr}.txt",
        "control": f"data/perbase_error_by_chr/{wildcards.genome}/{control}_{wildcards.strand}/{control}_{wildcards.strand}_{wildcards.chr}.txt"
    }


rule calculate_reactivity:
    """Calculate reactivity by comparing treatment to control samples (per chromosome)."""
    input:
        unpack(get_reactivity_inputs)
    output:
        reactivity="data/reactivity/{genome}/{sample}_{strand}_{chr}.txt.gz"
    log:
        "logs/reactivity/{genome}/{sample}_{strand}_{chr}.log"
    wildcard_constraints:
        strand="for|rev"
    shell:
        """
        workflow/scripts/Calculate_reactivity.sh \
            -t {input.treatment} \
            -c {input.control} \
            -o {output.reactivity} \
            2>&1 | tee {log}
        """
