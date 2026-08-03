# Phase 2b: Per-Base Pore Model Signal Deviation
# Only active when pod5 input is available (Uncalled4 TSV exists).
# Produces files in the same 5-column format as perbase_error so downstream
# phases can reuse Calculate_reactivity.sh and all phase 4 scripts.
#
# Input: strand-specific TSV from uncalled4_convert_tsv (phase 0).
# Produced by converting the Uncalled4 BAM (no re-alignment); pre-filtered
# to one strand by samtools before uncalled4 convert ran.
#
# Note: split_signal_by_chr_{genome} rules are generated per-genome in
# genome_specific_rules.smk (by CONFIG.sh) to avoid output-block wildcard issues.


rule perbase_signal_deviation:
    """Compute per-base pore model signal deviation from Uncalled4 DTW TSV.
    Parses dtw.model_diff (observed - expected pore model current) per reference
    position and outputs the 5-column format (chr, pos, nt, cov, mean_deviation).
    """
    input:
        tsv="data/uncalled4_tsv/{genome}/{raw_sample}_{strand}.tsv",
        genome="resources/genomes/{genome}.fa"
    output:
        dev="data/perbase_signal/{genome}/{raw_sample}_{strand}.txt.gz"
    log:
        "logs/perbase_signal/{genome}/{raw_sample}_{strand}.log"
    benchmark:
        "benchmarks/phase2b/perbase_signal_deviation/{genome}/{raw_sample}_{strand}.tsv"
    wildcard_constraints:
        strand="for|rev"
    shell:
        """
        python3 workflow/scripts/perbase_signal_deviation.py \
            -i {input.tsv} \
            -g {input.genome} \
            -o {output.dev} \
            2>&1 | tee {log}
        """
