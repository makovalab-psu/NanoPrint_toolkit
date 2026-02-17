#!/bin/bash

# Calculate Mean Reactivity: Average reactivity in genomic windows
# Computes mean reactivity per window without significance filtering

set -euo pipefail

# Default values
WINDOW_SIZE=1000

# Usage function
usage() {
    cat << EOF
Usage: $(basename "$0") -i <reactivity.txt.gz> -o <output.bg> -g <genome.fa.fai> [-w window] [-T tmpdir]

Calculate mean reactivity in genomic windows.

Required arguments:
    -i    Input reactivity file (gzipped, from Calculate_reactivity.sh)
    -o    Output bedGraph file (mean reactivity per window)
    -g    Chromosome sizes file (FAI format from samtools faidx)

Optional arguments:
    -w    Window size in bp (default: 1000)
    -T    Temporary directory (default: same directory as output)
    -h    Show this help message

Input format (4 columns, tab-separated):
    1. Chromosome name
    2. Position (1-based)
    3. Nucleotide identity
    4. Reactivity

Output format: bedGraph (4 columns)
    1. Chromosome name
    2. Start (0-based)
    3. End (1-based)
    4. Mean reactivity in window

Example:
    $(basename "$0") -i reactivity_chr1.txt.gz -o mean_chr1.bg -g genome.fa.fai -w 1000
EOF
    exit 1
}

# Parse arguments
INPUT=""
OUTPUT=""
GENOME=""
TMP_DIR=""

while getopts "i:o:g:w:T:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        g) GENOME="$OPTARG" ;;
        w) WINDOW_SIZE="$OPTARG" ;;
        T) TMP_DIR="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Validate required arguments
if [[ -z "$INPUT" || -z "$OUTPUT" || -z "$GENOME" ]]; then
    echo "Error: Missing required arguments" >&2
    usage
fi

# Check input files exist
if [[ ! -f "$INPUT" ]]; then
    echo "Error: Input file not found: $INPUT" >&2
    exit 1
fi

if [[ ! -f "$GENOME" ]]; then
    echo "Error: Genome sizes file not found: $GENOME" >&2
    exit 1
fi

# Check bedtools is available
if ! command -v bedtools &> /dev/null; then
    echo "Error: bedtools not found. Please install bedtools." >&2
    exit 1
fi

# Get output directory
OUT_DIR=$(dirname "$OUTPUT")
if [[ -z "$OUT_DIR" || "$OUT_DIR" == "." ]]; then
    OUT_DIR="$(pwd)"
fi
mkdir -p "$OUT_DIR"

# Set temp directory (default: in output directory)
if [[ -z "$TMP_DIR" ]]; then
    TMP_DIR="${OUT_DIR}/tmp_mean_bg_$$"
fi
mkdir -p "$TMP_DIR"

# Cleanup function
cleanup() {
    if [[ -d "$TMP_DIR" ]]; then
        rm -rf "$TMP_DIR"
    fi
}
trap cleanup EXIT

# Define temp file paths
DATA_BG="${TMP_DIR}/data.bg"
WINDOWS_BED="${TMP_DIR}/windows.bed"
FILTERED_GENOME="${TMP_DIR}/filtered_genome.txt"

echo "Decompressing and filtering input..."
# Decompress, remove 999999/-999999 markers, convert to bedGraph (0-based start)
gunzip -c "$INPUT" | awk -F'\t' '$4 != 999999 && $4 != -999999 { print $1, $2-1, $2, $4 }' OFS='\t' > "$DATA_BG"

DATA_COUNT=$(wc -l < "$DATA_BG" | tr -d ' ')
echo "Valid positions: $DATA_COUNT"

if [[ "$DATA_COUNT" -eq 0 ]]; then
    echo "Warning: No valid positions in input" >&2
    # Create empty output
    > "$OUTPUT"
    exit 0
fi

echo "Extracting chromosomes from input..."
# Get unique chromosomes from input
INPUT_CHRS=$(cut -f1 "$DATA_BG" | sort -u)
CHR_COUNT=$(echo "$INPUT_CHRS" | wc -l | tr -d ' ')
echo "Chromosomes in input: $CHR_COUNT ($(echo $INPUT_CHRS | tr '\n' ' '))"

# Filter genome file to only include chromosomes from input
echo "Filtering genome file to matching chromosomes..."
awk -v chrs="$INPUT_CHRS" '
BEGIN {
    n = split(chrs, arr, "\n")
    for (i in arr) valid[arr[i]] = 1
}
$1 in valid { print $1, $2 }
' OFS='\t' "$GENOME" > "$FILTERED_GENOME"

if [[ ! -s "$FILTERED_GENOME" ]]; then
    echo "Error: No matching chromosomes found between input and genome file" >&2
    exit 1
fi

echo "Creating genomic windows (${WINDOW_SIZE}bp)..."
bedtools makewindows -g "$FILTERED_GENOME" -w "$WINDOW_SIZE" > "$WINDOWS_BED"

echo "Calculating mean reactivity per window..."
# Use bedtools map to compute mean reactivity per window
# -a = windows, -b = reactivity positions
# -c 4 = operate on column 4 (reactivity)
# -o mean = compute mean
bedtools map -a "$WINDOWS_BED" -b "$DATA_BG" -c 4 -o mean | \
    awk 'BEGIN{OFS="\t"} {
        # Replace "." with 0 for windows with no data
        if ($4 == ".") $4 = 0
        print
    }' > "$OUTPUT"

echo ""
echo "Done: $OUTPUT"
