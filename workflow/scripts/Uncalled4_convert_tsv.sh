#!/bin/bash

# Uncalled4_convert_tsv.sh - Convert an Uncalled4 BAM to a strand-specific DTW TSV
# Pre-filters the BAM to one strand with samtools, then uses "uncalled4 convert"
# to extract DTW metrics without re-running signal alignment.
#
# Usage: Uncalled4_convert_tsv.sh -i <uncalled4.bam> -o <out.tsv> -s <for|rev> [-t <threads>]
#
# uncalled4 convert syntax (from official docs):
#   uncalled4 convert --bam-in <bam> --tsv-out <tsv> --tsv-cols "..."

set -euo pipefail

usage() {
    cat << EOF
Usage: $(basename "$0") -i <uncalled4.bam> -o <out.tsv> -s <for|rev> [-t <threads>]

Convert an Uncalled4 signal-alignment BAM to a strand-specific DTW TSV.
Pre-filters the BAM to the requested strand with samtools view, then calls
"uncalled4 convert" to extract DTW columns without re-running alignment.

Required arguments:
    -i    Input Uncalled4 BAM (data/uncalled4/{genome}/{sample}.bam)
    -o    Output TSV (dtw metrics per read × reference position)
    -s    Strand: for (forward, -F 0x10) or rev (reverse, -f 0x10)

Optional arguments:
    -t    Number of parallel processes for uncalled4 convert (default: 4)
    -h    Show this help message

TSV columns extracted from the BAM:
    dtw.current       Normalized mean read signal current (pA) per event
    dtw.current_sd    Signal current standard deviation
    dtw.start         Signal sample start index in raw trace
    dtw.length        Number of signal samples spanning this position
    dtw.model_diff    Model current - observed current (pA); positive = observed < expected
Note: dtw.base is NOT requested — it is not a valid layer in current uncalled4 versions
and causes "Invalid layer" ValueError. perbase_signal_deviation.py falls back to pysam
FASTA lookup for nucleotide identity when dtw.base is absent.

Example:
    $(basename "$0") \\
        -i data/uncalled4/genome/Sample01.bam \\
        -o data/uncalled4_tsv/genome/Sample01_for.tsv \\
        -s for \\
        -t 4
EOF
    exit 1
}

INPUT_BAM=""
OUTPUT_TSV=""
STRAND=""
THREADS=4

while getopts "i:o:s:t:h" opt; do
    case $opt in
        i) INPUT_BAM="$OPTARG" ;;
        o) OUTPUT_TSV="$OPTARG" ;;
        s) STRAND="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$INPUT_BAM" || -z "$OUTPUT_TSV" || -z "$STRAND" ]]; then
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

if ! command -v uncalled4 &> /dev/null; then
    echo "Error: uncalled4 not found in PATH. Install with: pip install uncalled4" >&2
    exit 1
fi

OUT_DIR=$(dirname "$OUTPUT_TSV")
mkdir -p "$OUT_DIR"

TMP_DIR="${OUT_DIR}/tmp_unc4conv_$$"
mkdir -p "$TMP_DIR"

cleanup() {
    [[ -d "$TMP_DIR" ]] && rm -rf "$TMP_DIR"
}
trap cleanup EXIT

# Forward: exclude reverse-strand reads (-F 0x10)
# Reverse: include only reverse-strand reads (-f 0x10)
if [[ "$STRAND" == "for" ]]; then
    STRAND_FLAG="-F 0x10"
    STRAND_LABEL="forward"
else
    STRAND_FLAG="-f 0x10"
    STRAND_LABEL="reverse"
fi

echo "=== Uncalled4 BAM → TSV Conversion ==="
echo "Input BAM:   $INPUT_BAM"
echo "Strand:      $STRAND_LABEL ($STRAND_FLAG)"
echo "Output TSV:  $OUTPUT_TSV"
echo "Processes:   $THREADS"
echo ""

# Step 1: Strand-filter the uncalled4 BAM
TMP_BAM="${TMP_DIR}/strand_filtered.bam"
echo "Pre-filtering BAM to $STRAND_LABEL strand..."
# shellcheck disable=SC2086
samtools view -b -h $STRAND_FLAG "$INPUT_BAM" > "$TMP_BAM"
samtools index "$TMP_BAM"
echo "Reads in strand-filtered BAM: $(samtools view -c "$TMP_BAM")"
echo ""

# Step 2: Convert BAM to TSV (no signal re-alignment needed)
# uncalled4 convert can return non-zero on partial failures; use || true + non-empty check.
echo "Converting to TSV..."
uncalled4 convert \
    --bam-in  "$TMP_BAM" \
    -p        "$THREADS" \
    --tsv-out "$OUTPUT_TSV" \
    --tsv-cols "dtw.current,dtw.current_sd,dtw.start,dtw.length,dtw.model_diff" || true

if [[ ! -f "$OUTPUT_TSV" ]]; then
    echo "Error: uncalled4 convert produced no TSV output" >&2
    exit 1
fi

LINE_COUNT=$(wc -l < "$OUTPUT_TSV")
if [[ "$LINE_COUNT" -le 1 ]]; then
    echo "Warning: output TSV is empty (header only or no output, $LINE_COUNT line(s))"
    echo "This is expected when the sample has few/no reads on this strand for this genome."
fi

echo ""
echo "Done. Output TSV: $OUTPUT_TSV"
echo "Row count (including header): $LINE_COUNT"
