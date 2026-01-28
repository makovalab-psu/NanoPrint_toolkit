#!/bin/bash

# Average Feature Annotation: Average annotations by distance, sample, and strand
# Chunks data for memory efficiency with large files

set -euo pipefail

# Usage function
usage() {
    cat << EOF
Usage: $(basename "$0") -i <input.txt.gz> -o <output.txt.gz> [-T tmpdir]

Average feature annotations by distance, sample (Treatment/Control), and strand.
Uses chunking for memory-efficient processing of large files.

Required arguments:
    -i    Input merged annotation file (gzipped)
    -o    Output averaged file (gzipped)

Optional arguments:
    -T    Temporary directory (default: output directory)
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

# Create output directory if needed
OUT_DIR=$(dirname "$OUTPUT")
if [[ -z "$OUT_DIR" || "$OUT_DIR" == "." ]]; then
    OUT_DIR="$(pwd)"
fi
mkdir -p "$OUT_DIR"

# Set temp directory
if [[ -z "$TMP_DIR" ]]; then
    TMP_DIR="${OUT_DIR}/tmp_average_$$"
fi
mkdir -p "$TMP_DIR"

CHUNK_DIR="$TMP_DIR/chunks"
mkdir -p "$CHUNK_DIR"

# Cleanup function
cleanup() {
    if [[ -d "$TMP_DIR" ]]; then
        rm -rf "$TMP_DIR"
    fi
}
trap cleanup EXIT INT TERM

echo "=== Average Feature Annotations ==="
echo "Input: $INPUT"
echo "Output: $OUTPUT"
echo ""

# Step 1: Chunk data by Distance, Sample, Strand
echo "Chunking data by distance, sample, and strand..."
gunzip -c "$INPUT" | awk -v chunk_dir="$CHUNK_DIR" '
BEGIN { FS = "\t"; OFS = "\t" }
NR == 1 { next }  # Skip header
{
    distance = $1
    sample = $5
    strand = $6

    # Create filename: distance_sample_strand.txt
    filename = chunk_dir "/" distance "_" sample "_" strand ".txt"
    print $0 >> filename
}
'

echo "Chunking complete"

# Step 2: Average each chunk
echo "Averaging chunks..."
AVERAGED="$TMP_DIR/averaged.txt"

# Write header
echo -e "Distance\tCoverage\tPerbase_error\tReactivity\tSample\tStrand" > "$AVERAGED"

for chunk_file in "$CHUNK_DIR"/*; do
    [[ -f "$chunk_file" ]] || continue

    # Extract metadata from filename: distance_sample_strand.txt
    basename=$(basename "$chunk_file" .txt)
    distance=$(echo "$basename" | cut -d'_' -f1)
    sample=$(echo "$basename" | cut -d'_' -f2)
    strand=$(echo "$basename" | cut -d'_' -f3)

    # Calculate averages
    awk -v dist="$distance" -v samp="$sample" -v str="$strand" '
    BEGIN { FS = "\t"; OFS = "\t" }
    {
        sum_cov += $2
        sum_err += $3
        count++

        # Only sum reactivity if present (Treatment rows)
        if ($4 != "") {
            sum_react += $4
            react_count++
        }
    }
    END {
        avg_cov = (count > 0) ? sum_cov / count : 0
        avg_err = (count > 0) ? sum_err / count : 0

        if (react_count > 0) {
            avg_react = sprintf("%.6f", sum_react / react_count)
        } else {
            avg_react = ""
        }

        printf "%s\t%.2f\t%.6f\t%s\t%s\t%s\n", dist, avg_cov, avg_err, avg_react, samp, str
    }
    ' "$chunk_file" >> "$AVERAGED"
done

# Step 3: Sort and compress output
echo "Sorting and compressing output..."
(head -n 1 "$AVERAGED" && tail -n +2 "$AVERAGED" | sort -t$'\t' -k1,1n -k5,5 -k6,6) | gzip -c > "$OUTPUT"

echo ""
echo "Done. Output: $OUTPUT"
