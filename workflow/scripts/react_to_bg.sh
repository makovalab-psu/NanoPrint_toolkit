#!/bin/bash

# Convert Reactivity file to BedGraph format
# Includes significance-based coloring

set -euo pipefail

# Usage function
usage() {
    cat << EOF
Usage: $(basename "$0") -i <reactivity.txt.gz> -o <output.bg> [-T tmpdir]

Convert reactivity data to bedGraph format with significance-based coloring.

Required arguments:
    -i    Input reactivity file (gzipped)
    -o    Output bedGraph file

Optional arguments:
    -T    Temporary directory (default: same directory as output)
    -h    Show this help message

Input format (4 columns, tab-separated):
    1. Chromosome name
    2. Position (1-based)
    3. Nucleotide identity
    4. Reactivity

Output format: UCSC bedGraph (4 columns: chr, start, end, reactivity)
    Header contains significance thresholds for reference:
    Thresholds: ns, p<0.05, p<0.01, p<0.001, p<1e-04
    Colors:     grey, black, #FF8C00, red, #810000

Example:
    $(basename "$0") -i reactivity_chr1.txt.gz -o reactivity_chr1.bg
EOF
    exit 1
}

# Parse arguments
INPUT=""
OUTPUT=""
TMP_DIR=""

while getopts "i:o:T:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        T) TMP_DIR="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Validate required arguments
if [[ -z "$INPUT" || -z "$OUTPUT" ]]; then
    echo "Error: Missing required arguments" >&2
    usage
fi

# Check input file exists
if [[ ! -f "$INPUT" ]]; then
    echo "Error: Input file not found: $INPUT" >&2
    exit 1
fi


# Get output directory (use current dir if output has no path)
OUT_DIR=$(dirname "$OUTPUT")
if [[ -z "$OUT_DIR" || "$OUT_DIR" == "." ]]; then
    OUT_DIR="$(pwd)"
fi
mkdir -p "$OUT_DIR"

# Set temp directory (default: in output directory)
if [[ -z "$TMP_DIR" ]]; then
    TMP_DIR="${OUT_DIR}/tmp_bg_$$"
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
FILTERED="${TMP_DIR}/filtered.txt"
ABS_VALUES="${TMP_DIR}/abs_values.txt"
SAMPLED="${TMP_DIR}/sampled.txt"
THRESHOLDS="${TMP_DIR}/thresholds.txt"

echo "Decompressing and filtering input..."
# Decompress and remove 999999/-999999 rows
gunzip -c "$INPUT" | awk -F'\t' '$4 != 999999 && $4 != -999999' > "$FILTERED"

echo "Extracting null distribution (absolute value of negative reactivity values)..."
# Negative reactivity values represent random chance - use as null distribution
# Take absolute value of negative values only
awk -F'\t' '$4 < 0 { print -$4 }' "$FILTERED" > "$ABS_VALUES"

NULL_COUNT=$(wc -l < "$ABS_VALUES" | tr -d ' ')
echo "Null distribution size: $NULL_COUNT negative values"

# If no negative values exist, skip threshold computation and write all data
if [[ "$NULL_COUNT" -eq 0 ]]; then
    echo "No negative reactivity values found - skipping significance threshold computation"
    echo "Writing bedGraph output (all data, no significance thresholds)..."
    {
        echo "# Reactivity bedGraph - No significance thresholds (null distribution empty)"
        echo ""
        awk -F'\t' '{ print $1, $2-1, $2, $4 }' OFS='\t' "$FILTERED"
    } > "$OUTPUT"
    echo "Done: $OUTPUT"
    exit 0
fi

echo "Downsampling for threshold calculation..."
# Randomly downsample to 1,000,000 rows (or all if fewer)
if [[ "$NULL_COUNT" -gt 1000000 ]]; then
    shuf -n 1000000 "$ABS_VALUES" | sort -n > "$SAMPLED"
    SAMPLE_SIZE=1000000
else
    sort -n "$ABS_VALUES" > "$SAMPLED"
    SAMPLE_SIZE="$NULL_COUNT"
fi
echo "Sample size: $SAMPLE_SIZE"

echo "Calculating significance thresholds from null distribution..."
# Thresholds based on percentiles of null distribution (negative values)
# p < 0.05: value exceeds 95% of null = 95th percentile
# p < 0.01: value exceeds 99% of null = 99th percentile
# p < 0.001: value exceeds 99.9% of null = 99.9th percentile
# p < 1e-04: value exceeds 99.99% of null = 99.99th percentile

# Use R's quantile function for accurate percentile calculation
# Using system R (/usr/local/bin/Rscript) to avoid conda library conflicts
read THRESH_05 THRESH_01 THRESH_001 THRESH_0001 <<< $(Rscript --vanilla -e "
x <- scan('$SAMPLED', quiet=TRUE)
q <- quantile(x, probs=c(0.95, 0.99, 0.999, 0.9999), na.rm = TRUE)
cat(q, sep=' ')
")

echo ""
echo "=== Significance Thresholds (from null distribution) ==="
echo "  ns (grey):           reactivity < $THRESH_05"
echo "  p < 0.05 (black):    reactivity >= $THRESH_05"
echo "  p < 0.01 (#FF8C00):  reactivity >= $THRESH_01"
echo "  p < 0.001 (red):     reactivity >= $THRESH_001"
echo "  p < 1e-04 (#810000): reactivity >= $THRESH_0001"
echo ""

echo "Writing bedGraph output..."
# Write bedGraph with multiple tracks for each significance level
{
    # Header comments with threshold information
    echo "# Reactivity bedGraph - Significance thresholds (from null distribution)"
    echo "# ns (grey):           reactivity < $THRESH_05"
    echo "# p < 0.05 (black):    reactivity >= $THRESH_05"
    echo "# p < 0.01 (#FF8C00):  reactivity >= $THRESH_01"
    echo "# p < 0.001 (red):     reactivity >= $THRESH_001"
    echo "# p < 1e-04 (#810000): reactivity >= $THRESH_0001"
    echo ""
    
    # Output all data with bedGraph coordinates (0-based start, 1-based end)
    awk -F'\t' '{ print $1, $2-1, $2, $4 }' OFS='\t' "$FILTERED"

} > "$OUTPUT"

echo "Done: $OUTPUT"
