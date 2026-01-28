#!/bin/bash

# Calculate React Density: Count permanganate reactive nucleotides in genomic windows
# Filters by significance threshold and counts per window

set -euo pipefail

# Default values
WINDOW_SIZE=1000
SIG_LEVEL=4

# Usage function
usage() {
    cat << EOF
Usage: $(basename "$0") -i <input.bg> -o <output.bg> -g <chrom.sizes.fai> [-w window] [-p sig_level] [-T tmpdir]

Calculate density of permanganate reactive nucleotides in genomic windows.

Required arguments:
    -i    Input bedGraph file (from react_to_bg.sh, with significance header)
    -o    Output bedGraph file (counts per window)
    -g    Chromosome sizes file (FAI format from samtools faidx)

Optional arguments:
    -w    Window size in bp (default: 1000)
    -p    Significance level filter (default: 4)
              1 = p <= 0.05
              2 = p <= 0.01
              3 = p <= 0.001
              4 = p <= 0.0001
    -T    Temporary directory (default: same directory as output)
    -h    Show this help message

Input format: bedGraph with significance thresholds in header
    Header example:
    # p < 0.05 (black):    reactivity >= 0.002886645
    # p < 0.01 (#FF8C00):  reactivity >= 0.004327449
    ...

Output format: bedGraph (5 columns)
    1. Chromosome name
    2. Start (0-based)
    3. End (1-based)
    4. Count of reactive nucleotides in window
    5. Sum of significant reactivity signal in window (after filtering)

Example:
    $(basename "$0") -i react_chr1.bg -o density_chr1.bg -g genome.fa.fai -w 1000 -p 4
EOF
    exit 1
}

# Parse arguments
INPUT=""
OUTPUT=""
GENOME=""
TMP_DIR=""

while getopts "i:o:g:w:p:T:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        g) GENOME="$OPTARG" ;;
        w) WINDOW_SIZE="$OPTARG" ;;
        p) SIG_LEVEL="$OPTARG" ;;
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

# Validate significance level
if [[ ! "$SIG_LEVEL" =~ ^[1-4]$ ]]; then
    echo "Error: Significance level (-p) must be 1, 2, 3, or 4" >&2
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
    TMP_DIR="${OUT_DIR}/tmp_dens_$$"
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
FILTERED_BG="${TMP_DIR}/filtered.bg"
WINDOWS_BED="${TMP_DIR}/windows.bed"
FILTERED_GENOME="${TMP_DIR}/filtered_genome.txt"

# Map significance level to header pattern
case $SIG_LEVEL in
    1) SIG_PATTERN="p < 0.05" ;;
    2) SIG_PATTERN="p < 0.01" ;;
    3) SIG_PATTERN="p < 0.001" ;;
    4) SIG_PATTERN="p < 1e-04" ;;
esac

echo "Parsing significance threshold from header..."
# Extract threshold value from header
# Header format: # p < 0.05 (black):    reactivity >= 0.002886645
THRESHOLD=$(grep "$SIG_PATTERN" "$INPUT" | grep -oE '[0-9]+\.[0-9e+-]+$' | head -1)

if [[ -z "$THRESHOLD" ]]; then
    echo "Error: Could not find threshold for '$SIG_PATTERN' in header" >&2
    echo "Header contents:" >&2
    head -10 "$INPUT" >&2
    exit 1
fi

echo "Significance level: $SIG_PATTERN"
echo "Reactivity threshold: >= $THRESHOLD"

echo "Filtering bedGraph by threshold..."
# Filter out header and values below threshold
awk -v thresh="$THRESHOLD" '
    !/^#/ && !/^$/ && $4 >= thresh { print }
' "$INPUT" > "$FILTERED_BG"

FILTERED_COUNT=$(wc -l < "$FILTERED_BG" | tr -d ' ')
echo "Filtered positions: $FILTERED_COUNT"

if [[ "$FILTERED_COUNT" -eq 0 ]]; then
    echo "Warning: No positions passed the significance threshold" >&2
fi

echo "Extracting chromosomes from input..."
# Get unique chromosomes from input bedGraph (skip header lines)
INPUT_CHRS=$(grep -v "^#" "$INPUT" | grep -v "^$" | cut -f1 | sort -u)
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
    echo "Input chromosomes: $INPUT_CHRS" >&2
    echo "Genome file chromosomes:" >&2
    cut -f1 "$GENOME" | head -10 >&2
    exit 1
fi

echo "Creating genomic windows (${WINDOW_SIZE}bp)..."
# Create windows file using bedtools makewindows
bedtools makewindows -g "$FILTERED_GENOME" -w "$WINDOW_SIZE" > "$WINDOWS_BED"

echo "Calculating coverage and sum per window..."
# Use bedtools map to count reactive nucleotides and sum reactivity per window
# -a = windows, -b = filtered reactive positions
# -c 4 = operate on column 4 (reactivity)
# -o count,sum = output count and sum
bedtools map -a "$WINDOWS_BED" -b "$FILTERED_BG" -c 4 -o count,sum | \
    awk 'BEGIN{OFS="\t"} {
        # Replace "." with 0 for windows with no overlapping features
        if ($4 == ".") $4 = 0
        if ($5 == ".") $5 = 0
        print
    }' > "$OUTPUT"


echo ""
echo "Done: $OUTPUT"
