# Phase 2b: Per-Base Pore Model Signal Deviation
# Only active when pod5 input is available (Uncalled4 TSV exists).
# Produces a 10-column file: columns 1-5 match perbase_error (chr, pos, nt, cov,
# mean), followed by four quantiles and the mean squared deviation. Downstream
# phases still reuse Calculate_reactivity.sh and all phase 4 scripts unchanged,
# because those select the value column with -f; phase 3b passes -f 10.
#
# Input: the Uncalled4 BAM itself. The DTW tags are read straight out of it with
# uncalled4's own decoder, so uncalled4_convert_tsv is no longer on this path -
# it stays a rule for small datasets and for producing test oracles. This removes
# the ~570 GB per-sample-strand TSV and the strand-filtered BAM copy it needed;
# the strand is now selected in memory, from the -s flag.
#
# Note: split_signal_by_chr_{genome} rules are generated per-genome in
# genome_specific_rules.smk (by CONFIG.sh) to avoid output-block wildcard issues.


rule perbase_signal_deviation:
    """Compute per-base pore model signal deviation from the Uncalled4 BAM.
    Parses dtw.model_diff (observed - expected pore model current) per reference
    position and writes 10 columns: chr, pos, nt, cov, mean_deviation, then the
    Q25/Q75/Q025/Q975 quantiles and the mean squared deviation (column 10, the
    column phase 3b selects with -f 10).
    """
    input:
        bam=uncalled4_bam,
        bai=uncalled4_bai,
        genome="resources/genomes/{genome}.fa"
    output:
        dev="data/perbase_signal/{genome}/{raw_sample}_{strand}.txt.gz"
    params:
        # Sets peak memory: one window's observations are held at once, so roughly
        # window x per-strand coverage x 12 bytes (~120 MB at 10 kb and 1000x).
        # Lower it for deeper data; it cannot change the result, only the peak.
        window=config.get("signal_window_size", 10000)
    log:
        "logs/perbase_signal/{genome}/{raw_sample}_{strand}.log"
    benchmark:
        "benchmarks/phase2b/perbase_signal_deviation/{genome}/{raw_sample}_{strand}.tsv"
    wildcard_constraints:
        strand="for|rev"
    shell:
        """
        python3 workflow/scripts/perbase_signal_deviation.py \
            -i {input.bam} \
            -g {input.genome} \
            -s {wildcards.strand} \
            -w {params.window} \
            -o {output.dev} \
            2>&1 | tee {log}
        """
