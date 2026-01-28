# Phase 3: Reactivity Calculation
# Rules for calculating reactivity (treatment - control) per chromosome


rule calculate_reactivity:
    """Calculate reactivity by comparing treatment to control samples (per chromosome)."""
    input:
        treatment="data/perbase_error_by_chr/{genome}/{treatment_sample}_{strand}/{treatment_sample}_{strand}_{chr}.txt",
        control="data/perbase_error_by_chr/{genome}/{control_sample}_{strand}/{control_sample}_{strand}_{chr}.txt"
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
