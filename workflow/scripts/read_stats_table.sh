#!/bin/bash

# Read Stats Table: Combine per-sample read statistics into a single CSV table

set -euo pipefail

usage() {
    cat << EOF
Usage: $(basename "$0") -o <output.csv> <input1.txt> [input2.txt ...]

Combine per-sample read statistics into a single CSV table.

Required arguments:
    -o    Output CSV file

Positional arguments:
    One or more read_stats txt files (tab-delimited, output of Read_stats.sh)

Output format (CSV):
    Sample, Giga_bp, Reads_million, N50, Q50

Example:
    $(basename "$0") -o tables/read_stats_table.csv tables/read_stats/sample1.txt tables/read_stats/sample2.txt
EOF
    exit 1
}

OUTPUT=""

while getopts "o:h" opt; do
    case $opt in
        o) OUTPUT="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done
shift $((OPTIND - 1))

INPUTS=("$@")

if [[ -z "$OUTPUT" || ${#INPUTS[@]} -eq 0 ]]; then
    echo "Error: Missing required arguments" >&2
    usage
fi

mkdir -p "$(dirname "$OUTPUT")"

echo "Combining ${#INPUTS[@]} read stats files..."

# Write CSV header
echo "Sample,Giga_bp,Reads_million,N50,Q50" > "$OUTPUT"

# Process each input file: skip header, convert tab-delimited to CSV
for f in "${INPUTS[@]}"; do
    if [[ ! -f "$f" ]]; then
        echo "Warning: File not found: $f" >&2
        continue
    fi
    tail -n +2 "$f" | awk 'BEGIN { OFS = "," } { print $1, $2, $3, $4, $5 }' >> "$OUTPUT"
done

echo "Done: $OUTPUT"
