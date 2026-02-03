#!/bin/bash

# Average Feature Annotation: Average annotations by distance, sample, and strand
# Uses single-pass awk with associative arrays for efficiency

set -euo pipefail

# Usage function
usage() {
    cat << EOF
Usage: $(basename "$0") -i <input.txt.gz> -o <output.txt.gz>

Average feature annotations by distance, sample (Treatment/Control), and strand.
Uses single-pass processing with associative arrays for efficiency.

Required arguments:
    -i    Input merged annotation file (gzipped)
    -o    Output averaged file (gzipped)

Optional arguments:
    -h    Show this help message

Input format (tab-delimited, gzipped):
    Distance, Coverage, Perbase_error, Reactivity, Sample, Strand

Output format (tab-delimited, gzipped):
    Distance, Coverage, Perbase_error, Reactivity, Sample, Strand
    (averaged across all features for each Distance/Sample/Strand combination)

Example:
    $(basename "$0") -i merged_annotations.txt.gz -o averaged_annotations.txt.gz
EOF
    exit 1
}

# Parse arguments
INPUT=""
OUTPUT=""

while getopts "i:o:T:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        T) : ;;  # Ignore -T for backwards compatibility
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

# Create output directory if needed
OUT_DIR=$(dirname "$OUTPUT")
if [[ -z "$OUT_DIR" || "$OUT_DIR" == "." ]]; then
    OUT_DIR="$(pwd)"
fi
mkdir -p "$OUT_DIR"

echo "=== Average Feature Annotations ==="
echo "Input: $INPUT"
echo "Output: $OUTPUT"
echo ""

# Single-pass averaging using awk associative arrays
# Keys are "distance\tsample\tstrand", values are accumulated sums and counts
echo "Averaging annotations (single-pass)..."

gunzip -c "$INPUT" | awk '
BEGIN {
    FS = "\t"
    OFS = "\t"
}
NR == 1 { next }  # Skip header
{
    # Build key from distance, sample, strand
    key = $1 SUBSEP $5 SUBSEP $6

    sum_cov[key] += $2
    sum_err[key] += $3
    count[key]++

    # Only accumulate reactivity if present (Treatment rows)
    if ($4 != "") {
        sum_react[key] += $4
        react_count[key]++
    }
}
END {
    # Output header
    print "Distance\tCoverage\tPerbase_error\tReactivity\tSample\tStrand"

    for (key in count) {
        # Split key back into components
        split(key, k, SUBSEP)
        dist = k[1]
        samp = k[2]
        str = k[3]

        avg_cov = sum_cov[key] / count[key]
        avg_err = sum_err[key] / count[key]

        if (react_count[key] > 0) {
            avg_react = sprintf("%.6f", sum_react[key] / react_count[key])
        } else {
            avg_react = ""
        }

        printf "%s\t%.2f\t%.6f\t%s\t%s\t%s\n", dist, avg_cov, avg_err, avg_react, samp, str
    }
}
' | {
    # Read header first, then sort the rest numerically by distance
    IFS= read -r header
    echo "$header"
    sort -t$'\t' -k1,1n -k5,5 -k6,6
} | gzip -c > "$OUTPUT"

echo ""
echo "Done. Output: $OUTPUT"
