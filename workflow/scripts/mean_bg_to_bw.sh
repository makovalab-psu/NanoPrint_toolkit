#!/bin/bash

# Mean BedGraph to BigWig: Convert a mean reactivity bedGraph to bigWig format
# Simple conversion without significance header parsing

set -euo pipefail

# Usage function
usage() {
    cat << EOF
Usage: $(basename "$0") -i <input.bg> -o <output.bw> -g <genome.fa.fai> [-T tmpdir]

Convert a mean reactivity bedGraph file to bigWig format.

Required arguments:
    -i    Input bedGraph file (4-column: chr, start, end, mean_reactivity)
    -o    Output bigWig file
    -g    Chromosome sizes file (FAI format from samtools faidx)

Optional arguments:
    -T    Temporary directory (default: same directory as output)
    -h    Show this help message

Dependencies:
    - bedGraphToBigWig (UCSC tools)

Example:
    $(basename "$0") -i mean_merged.bg -o mean_merged.bw -g genome.fa.fai
EOF
    exit 1
}

# Parse arguments
INPUT=""
OUTPUT=""
GENOME_FAI=""
TMP_DIR=""

while getopts "i:o:g:T:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        g) GENOME_FAI="$OPTARG" ;;
        T) TMP_DIR="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Validate required arguments
if [[ -z "$INPUT" || -z "$OUTPUT" || -z "$GENOME_FAI" ]]; then
    echo "Error: Missing required arguments" >&2
    usage
fi

# Check input files exist
if [[ ! -f "$INPUT" ]]; then
    echo "Error: Input file not found: $INPUT" >&2
    exit 1
fi

if [[ ! -f "$GENOME_FAI" ]]; then
    echo "Error: Genome sizes file not found: $GENOME_FAI" >&2
    exit 1
fi

# Check required tools
if ! command -v bedGraphToBigWig &> /dev/null; then
    echo "Error: bedGraphToBigWig not found. Please install UCSC tools." >&2
    exit 1
fi

# Create output directory if needed
OUT_DIR=$(dirname "$OUTPUT")
if [[ -z "$OUT_DIR" || "$OUT_DIR" == "." ]]; then
    OUT_DIR="$(pwd)"
fi
mkdir -p "$OUT_DIR"

# Set temp directory (default: in output directory)
if [[ -z "$TMP_DIR" ]]; then
    TMP_DIR="${OUT_DIR}/tmp_mean_bw_$$"
fi
mkdir -p "$TMP_DIR"

# Cleanup function
cleanup() {
    if [[ -d "$TMP_DIR" ]]; then
        rm -rf "$TMP_DIR"
    fi
}
trap cleanup EXIT

# Create chrom.sizes from fai (columns 1 and 2)
CHROM_SIZES="${TMP_DIR}/chrom.sizes"
cut -f1,2 "$GENOME_FAI" > "$CHROM_SIZES"

# Sort bedGraph (required by bedGraphToBigWig)
echo "Sorting bedGraph..."
SORTED_BG="${TMP_DIR}/sorted.bg"
sort -k1,1 -k2,2n "$INPUT" > "$SORTED_BG"

# Convert to bigWig
echo "Converting to bigWig..."
bedGraphToBigWig "$SORTED_BG" "$CHROM_SIZES" "$OUTPUT"

echo "Done: $OUTPUT"
