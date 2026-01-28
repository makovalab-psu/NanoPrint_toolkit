#!/bin/bash

# Perbase Error: Calculate per-base error rates from BAM alignments
# Uses samtools mpileup and custom AWK script to compute error probabilities

set -euo pipefail

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Usage function
usage() {
    cat << EOF
Usage: $(basename "$0") <--rev|--for> -i <input.bam> -g <genome.fasta> -o <output.txt.gz> [-T tmpdir]

Calculate per-base error rates from filtered BAM alignments.

Positional argument (required, must be first):
    --rev    Process reverse strand reads only (FLAG 0x10 set)
    --for    Process forward strand reads only (FLAG 0x10 not set)

Required arguments:
    -i    Input BAM file (filtered/merged alignments)
    -g    Reference genome FASTA file
    -o    Output file (gzipped, *.txt.gz)

Optional arguments:
    -T    Temporary directory (default: same directory as output)
    -h    Show this help message

Output format (tab-delimited, gzipped):
    Column 1: Chromosome name
    Column 2: Position in chromosome (1-based)
    Column 3: Nucleotide identity
    Column 4: Coverage
    Column 5: Per-base error probability

Example:
    $(basename "$0") --for -i sample_filtered.bam -g reference.fasta -o sample_forward.txt.gz
    $(basename "$0") --rev -i sample_filtered.bam -g reference.fasta -o sample_reverse.txt.gz
EOF
    exit 1
}

# Check for strand argument first (positional)
if [[ $# -lt 1 ]]; then
    echo "Error: Missing strand argument (--rev or --for)" >&2
    usage
fi

STRAND_ARG="$1"
shift

if [[ "$STRAND_ARG" == "--rev" ]]; then
    SAMTOOLS_VIEW_FLAGS="-f 0x10"
    STRAND_DESC="reverse"
elif [[ "$STRAND_ARG" == "--for" ]]; then
    SAMTOOLS_VIEW_FLAGS="-F 0x10"
    STRAND_DESC="forward"
elif [[ "$STRAND_ARG" == "-h" ]]; then
    usage
else
    echo "Error: First argument must be --rev or --for (got: $STRAND_ARG)" >&2
    usage
fi

# Parse remaining arguments
INPUT_BAM=""
INPUT_FASTA=""
OUTPUT=""
TMP_DIR=""

while getopts "i:g:o:T:h" opt; do
    case $opt in
        i) INPUT_BAM="$OPTARG" ;;
        g) INPUT_FASTA="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        T) TMP_DIR="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Validate required arguments
if [[ -z "$INPUT_BAM" || -z "$INPUT_FASTA" || -z "$OUTPUT" ]]; then
    echo "Error: Missing required arguments" >&2
    usage
fi

# Check input files exist
if [[ ! -f "$INPUT_BAM" ]]; then
    echo "Error: Input BAM file not found: $INPUT_BAM" >&2
    exit 1
fi

if [[ ! -f "$INPUT_FASTA" ]]; then
    echo "Error: Reference genome file not found: $INPUT_FASTA" >&2
    exit 1
fi

# Check AWK script exists
AWK_SCRIPT="${SCRIPT_DIR}/perbase_error.awk"
if [[ ! -f "$AWK_SCRIPT" ]]; then
    echo "Error: AWK script not found: $AWK_SCRIPT" >&2
    exit 1
fi

# Check required tools are available
if ! command -v samtools &> /dev/null; then
    echo "Error: samtools not found. Please install samtools." >&2
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
    TMP_DIR="${OUT_DIR}/tmp_perbase_$$"
fi
mkdir -p "$TMP_DIR"

# Define temporary BAM file
TMP_BAM="${TMP_DIR}/strand_filtered.bam"

# Cleanup function
cleanup() {
    if [[ -d "$TMP_DIR" ]]; then
        rm -rf "$TMP_DIR"
    fi
}
trap cleanup EXIT INT TERM

echo "=== Per-base Error Calculation ==="
echo "Input BAM: $INPUT_BAM"
echo "Reference: $INPUT_FASTA"
echo "Strand: $STRAND_DESC ($STRAND_ARG)"
echo "Output: $OUTPUT"
echo ""

# Step 1: Filter BAM by strand
echo "Filtering reads by strand ($STRAND_DESC)..."
samtools view -b -h $SAMTOOLS_VIEW_FLAGS "$INPUT_BAM" > "$TMP_BAM"

# Step 2: Run mpileup and calculate per-base error
echo "Calculating per-base error rates..."
samtools mpileup -f "$INPUT_FASTA" -B -Q 0 "$TMP_BAM" | awk -f "$AWK_SCRIPT" | gzip -c > "$OUTPUT"

echo ""
echo "Done. Output: $OUTPUT"
