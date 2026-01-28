#!/bin/bash

# Split by Chromosome: Split a tab-delimited file by chromosome
# Enables parallelization of downstream analysis

set -euo pipefail

# Usage function
usage() {
    cat << EOF
Usage: $(basename "$0") -i <input_file>

Split a tab-delimited file by chromosome (first column).
Output files are created in the same directory as input.

Required arguments:
    -i    Input file (tab-delimited with chromosome in column 1)
          Supported formats: .txt, .txt.gz, .bed, .bed.gz, .bg, .bg.gz

Optional arguments:
    -h    Show this help message

Output:
    Creates <prefix>_<chr>.<ext> for each chromosome in the same directory.
    Prefix is derived from input filename (without extension).

Example:
    $(basename "$0") -i data/perbase_error/sample_for.txt.gz
    # Creates: data/perbase_error/sample_for_chr1.txt, sample_for_chr2.txt, ...
EOF
    exit 1
}

# Parse arguments
INPUT=""

while getopts "i:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Validate required arguments
if [[ -z "$INPUT" ]]; then
    echo "Error: Missing required arguments" >&2
    usage
fi

# Check input file exists
if [[ ! -f "$INPUT" ]]; then
    echo "Error: Input file not found: $INPUT" >&2
    exit 1
fi

# Determine output prefix and extension by removing suffix
INPUT_DIR=$(dirname "$INPUT")
INPUT_BASE=$(basename "$INPUT")

COMPRESSED=false
if [[ "$INPUT_BASE" == *.gz ]]; then
    COMPRESSED=true
    INPUT_BASE="${INPUT_BASE%.gz}"
fi

if [[ "$INPUT_BASE" == *.txt ]]; then
    PREFIX="${INPUT_BASE%.txt}"
    EXT="txt"
elif [[ "$INPUT_BASE" == *.bed ]]; then
    PREFIX="${INPUT_BASE%.bed}"
    EXT="bed"
elif [[ "$INPUT_BASE" == *.bg ]]; then
    PREFIX="${INPUT_BASE%.bg}"
    EXT="bg"
else
    echo "Error: Input file must have .txt, .bed, or .bg extension (optionally gzipped)" >&2
    exit 1
fi

OUT_PREFIX="${INPUT_DIR}/${PREFIX}_"

echo "=== Split by Chromosome ==="
echo "Input: $INPUT"
echo "Output prefix: $OUT_PREFIX"
echo "Output extension: $EXT"
echo ""

# Split file by chromosome
echo "Splitting file by chromosome..."
if [[ "$COMPRESSED" == true ]]; then
    gunzip -c "$INPUT" | awk -v prefix="$OUT_PREFIX" -v ext="$EXT" '{
        output_file = prefix $1 "." ext
        print > output_file
    }'
else
    awk -v prefix="$OUT_PREFIX" -v ext="$EXT" '{
        output_file = prefix $1 "." ext
        print > output_file
    }' "$INPUT"
fi

echo ""
echo "Done."
