# Phase 3b: Signal Reactivity Calculation
# Computes treatment_signal_deviation - control_signal_deviation per chromosome.
# Reuses Calculate_reactivity.sh since the 5-column format is identical to
# perbase_error. Output paths use 'signal_reactivity' to avoid collisions.


def get_signal_reactivity_inputs(wildcards):
    """Get per-chromosome signal deviation files for treatment and control."""
    treatment = get_treatment(wildcards.sample)
    control = get_control(wildcards.sample)
    return {
        "treatment": (
            f"data/perbase_signal_by_chr/{wildcards.genome}"
            f"/{treatment}_{wildcards.strand}"
            f"/{treatment}_{wildcards.strand}_{wildcards.chr}.txt.gz"
        ),
        "control": (
            f"data/perbase_signal_by_chr/{wildcards.genome}"
            f"/{control}_{wildcards.strand}"
            f"/{control}_{wildcards.strand}_{wildcards.chr}.txt.gz"
        ),
    }


rule calculate_signal_reactivity:
    """Calculate signal reactivity: treatment_deviation - control_deviation (per chromosome)."""
    input:
        unpack(get_signal_reactivity_inputs)
    output:
        reactivity=wrap_output(
            "signal_reactivity",
            "data/signal_reactivity/{genome}/{sample}_{strand}_{chr}.txt.gz"
        )
    log:
        "logs/signal_reactivity/{genome}/{sample}_{strand}_{chr}.log"
    benchmark:
        "benchmarks/phase3b/calculate_signal_reactivity/{genome}/{sample}_{strand}_{chr}.tsv"
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
