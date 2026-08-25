#!/bin/bash

# Uncalled4_align_tsv.sh - Strand-specific Uncalled4 DTW alignment with TSV output
# Pre-filters the input BAM to one strand, then runs Uncalled4 with --tsv-out.
# The TSV contains dtw.model_diff (observed - model pore model current per base,
# in normalized units),
# which perbase_signal_deviation.py aggregates into a per-position signal deviation file.
#
# Usage: Uncalled4_align_tsv.sh -i <filtered.bam> -p <pod5_dir> -g <genome.fa> -o <out.tsv> -s <for|rev> [-t <threads>]
#
# Command syntax confirmed from js4004 Snakefile:
#   uncalled4 align --bam-in <bam> --ref <fa> --reads <pod5> --tsv-out <tsv> --tsv-cols "..."
# Parallelism via -p (official docs: default is 1 process).

set -euo pipefail

usage() {
    cat << EOF
Usage: $(basename "$0") -i <filtered.bam> -p <pod5_dir> -g <genome.fa> -o <out.tsv> -s <for|rev> [-t <threads>]

Run Uncalled4 DTW alignment on a strand-specific BAM subset and write TSV output.
The TSV will contain dtw.model_diff per reference position, used by
perbase_signal_deviation.py to compute per-base signal deviation.

Required arguments:
    -i    Input sequence-aligned BAM (filtered_alignments/{genome}/{sample}.bam)
    -p    Pod5 file directory (raw_data/{sample}/) or single pod5 file
    -g    Reference genome FASTA
    -o    Output TSV (dtw metrics per reference position)
    -s    Strand: for (forward, -F 0x10) or rev (reverse, -f 0x10)

Optional arguments:
    -t    Number of parallel processes for Uncalled4 (default: 4)
    -h    Show this help message

TSV columns requested from Uncalled4:
    dtw.current       Normalized mean read signal current (pA) per event
    dtw.current_sd    Signal current standard deviation
    dtw.start         Signal sample start index in raw trace
    dtw.length        Number of signal samples spanning this position
    dtw.model_diff    Observed - model current, normalized units (NOT pA);
                      positive = observed current higher than the pore model expects
    dtw.base          Binarized reference base (used as nucleotide; may be integer-encoded)

Example:
    $(basename "$0") \\
        -i data/filtered_alignments/genome/Sample01.bam \\
        -p raw_data/Sample01/ \\
        -g resources/genomes/genome.fa \\
        -o data/uncalled4_tsv/genome/Sample01_for.tsv \\
        -s for \\
        -t 4
EOF
    exit 1
}

INPUT_BAM=""
POD5_DIR=""
GENOME=""
OUTPUT_TSV=""
STRAND=""
THREADS=4

while getopts "i:p:g:o:s:t:h" opt; do
    case $opt in
        i) INPUT_BAM="$OPTARG" ;;
        p) POD5_DIR="$OPTARG" ;;
        g) GENOME="$OPTARG" ;;
        o) OUTPUT_TSV="$OPTARG" ;;
        s) STRAND="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$INPUT_BAM" || -z "$POD5_DIR" || -z "$GENOME" || -z "$OUTPUT_TSV" || -z "$STRAND" ]]; then
    echo "Error: Missing required arguments" >&2
    usage
fi

if [[ "$STRAND" != "for" && "$STRAND" != "rev" ]]; then
    echo "Error: -s must be 'for' or 'rev'" >&2
    usage
fi

if [[ ! -f "$INPUT_BAM" ]]; then
    echo "Error: Input BAM not found: $INPUT_BAM" >&2
    exit 1
fi

if [[ ! -e "$POD5_DIR" ]]; then
    echo "Error: Pod5 path not found: $POD5_DIR" >&2
    exit 1
fi

if [[ ! -f "$GENOME" ]]; then
    echo "Error: Genome FASTA not found: $GENOME" >&2
    exit 1
fi

if ! command -v uncalled4 &> /dev/null; then
    echo "Error: uncalled4 not found in PATH. Install with: pip install uncalled4" >&2
    exit 1
fi

OUT_DIR=$(dirname "$OUTPUT_TSV")
mkdir -p "$OUT_DIR"

TMP_DIR="${OUT_DIR}/tmp_unc4tsv_$$"
mkdir -p "$TMP_DIR"

cleanup() {
    [[ -d "$TMP_DIR" ]] && rm -rf "$TMP_DIR"
}
trap cleanup EXIT

# Forward: exclude reverse-strand reads (-F 0x10 = exclude reads where flag 0x10 is set)
# Reverse: include only reverse-strand reads (-f 0x10 = require flag 0x10)
if [[ "$STRAND" == "for" ]]; then
    STRAND_FLAG="-F 0x10"
    STRAND_LABEL="forward"
else
    STRAND_FLAG="-f 0x10"
    STRAND_LABEL="reverse"
fi

echo "=== Uncalled4 Signal Alignment (TSV) ==="
echo "Input BAM:   $INPUT_BAM"
echo "Pod5 path:   $POD5_DIR"
echo "Genome:      $GENOME"
echo "Strand:      $STRAND_LABEL ($STRAND_FLAG)"
echo "Output TSV:  $OUTPUT_TSV"
echo "Processes:   $THREADS"
echo ""

# Step 1: Pre-filter BAM to one strand
TMP_BAM="${TMP_DIR}/strand_filtered.bam"
echo "Pre-filtering BAM to $STRAND_LABEL strand..."
# shellcheck disable=SC2086
samtools view -b -h $STRAND_FLAG "$INPUT_BAM" > "$TMP_BAM"
samtools index "$TMP_BAM"
echo "Reads in strand-filtered BAM: $(samtools view -c "$TMP_BAM")"
echo ""

# Step 2: Run Uncalled4 with TSV output
# uncalled4 returns non-zero when any reads fail DTW (even if most succeed).
# Use || true and verify the output is non-empty (pattern from js4007/Snakefile).
echo "Running Uncalled4 DTW alignment..."
uncalled4 align \
    --bam-in  "$TMP_BAM" \
    --ref     "$GENOME" \
    --reads   "$POD5_DIR" \
    -p        "$THREADS" \
    --tsv-out "$OUTPUT_TSV" \
    --tsv-cols "dtw.current,dtw.current_sd,dtw.start,dtw.length,dtw.model_diff,dtw.base" || true

if [[ ! -f "$OUTPUT_TSV" ]]; then
    echo "Error: Uncalled4 produced no TSV output" >&2
    exit 1
fi

# Header line only = no data rows
LINE_COUNT=$(wc -l < "$OUTPUT_TSV")
if [[ "$LINE_COUNT" -le 1 ]]; then
    echo "Error: Uncalled4 TSV is empty (header only, $LINE_COUNT line(s))" >&2
    exit 1
fi

echo ""
echo "Done. Output TSV: $OUTPUT_TSV"
echo "Row count (including header): $LINE_COUNT"
