# Phase 3b: Signal Reactivity Calculation
# Computes treatment_mean_sq - control_mean_sq per chromosome, where mean_sq
# is mean(dtw.model_diff^2) = sum(dtw.model_diff^2)/N (pA^2, column 10 of
# perbase_signal files). Uses Calculate_reactivity.sh with -f 10 to select
# the mean squared deviation column from the 10-column perbase_signal format.
# Output paths use 'signal_reactivity' to avoid collisions with perbase_error reactivity.


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
            -f 10 \
            2>&1 | tee {log}
        """
