#!/bin/bash

# Merge Annotations: Merge chromosome-split annotation files
# Concatenates files preserving header from first file only

set -euo pipefail

# Usage function
usage() {
    cat << EOF
Usage: $(basename "$0") -o <output.txt.gz> <input1.txt.gz> <input2.txt.gz> ...

Merge chromosome-split annotation files into a single file.

Required arguments:
    -o    Output file (gzipped)

Positional arguments:
    Remaining arguments are input annotation files (gzipped)

Optional arguments:
    -h    Show this help message

Example:
    $(basename "$0") -o merged.txt.gz chr1.txt.gz chr2.txt.gz chr3.txt.gz
EOF
    exit 1
}

# Parse arguments
OUTPUT=""

while getopts "o:h" opt; do
    case $opt in
        o) OUTPUT="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Shift past the parsed options
shift $((OPTIND - 1))

# Remaining arguments are input files
INPUT_FILES=("$@")

# Validate required arguments
if [[ -z "$OUTPUT" ]]; then
    echo "Error: Missing required argument (-o)" >&2
    usage
fi

if [[ ${#INPUT_FILES[@]} -eq 0 ]]; then
    echo "Error: No input files specified" >&2
    usage
fi

# Create output directory if needed
OUT_DIR=$(dirname "$OUTPUT")
if [[ -n "$OUT_DIR" && "$OUT_DIR" != "." ]]; then
    mkdir -p "$OUT_DIR"
fi

echo "=== Merge Annotations ==="
echo "Output: $OUTPUT"
echo "Input files: ${#INPUT_FILES[@]}"
echo ""

# Merge files: header from first file, data from all files
{
    first=true
    for f in "${INPUT_FILES[@]}"; do
        if [[ ! -f "$f" ]]; then
            echo "Warning: File not found, skipping: $f" >&2
            continue
        fi

        if [[ "$first" == true ]]; then
            # First file: include header
            gunzip -c "$f"
            first=false
        else
            # Subsequent files: skip header
            gunzip -c "$f" | tail -n +2
        fi
    done
} | gzip -c > "$OUTPUT"

echo "Done. Output: $OUTPUT"
