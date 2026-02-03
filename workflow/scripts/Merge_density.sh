#!/bin/bash

# Merge Density: Merge chromosome-split density files in genome order
# Concatenates files in the order specified by the genome sizes file

set -euo pipefail

# Usage function
usage() {
    cat << EOF
Usage: $(basename "$0") -g <genome.fa.fai> -o <output.bg> <input1.bg> <input2.bg> ...

Merge chromosome-split density/bedgraph files in genome order.

Required arguments:
    -g    Genome sizes file (samtools faidx .fai format)
    -o    Output file

Positional arguments:
    Remaining arguments are input files to merge

Optional arguments:
    -h    Show this help message

The input files are merged in the order of chromosomes as they appear
in the genome sizes file (first column of .fai).

Example:
    $(basename "$0") -g genome.fa.fai -o merged.bg sample_chr1.bg sample_chr2.bg sample_chrX.bg
EOF
    exit 1
}

# Parse arguments
GENOME_FAI=""
OUTPUT=""

while getopts "g:o:h" opt; do
    case $opt in
        g) GENOME_FAI="$OPTARG" ;;
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
if [[ -z "$GENOME_FAI" || -z "$OUTPUT" ]]; then
    echo "Error: Missing required arguments (-g and -o)" >&2
    usage
fi

if [[ ${#INPUT_FILES[@]} -eq 0 ]]; then
    echo "Error: No input files specified" >&2
    usage
fi

# Check genome sizes file exists
if [[ ! -f "$GENOME_FAI" ]]; then
    echo "Error: Genome sizes file not found: $GENOME_FAI" >&2
    exit 1
fi

# Create output directory if needed
OUT_DIR=$(dirname "$OUTPUT")
if [[ -n "$OUT_DIR" && "$OUT_DIR" != "." ]]; then
    mkdir -p "$OUT_DIR"
fi

echo "=== Merge Density Files ==="
echo "Genome: $GENOME_FAI"
echo "Output: $OUTPUT"
echo "Input files: ${#INPUT_FILES[@]}"
echo ""

# Build parallel arrays of input files by chromosome (Bash 3.2 compatible)
declare -a CHR_NAMES=()
declare -a CHR_FILES=()
for f in "${INPUT_FILES[@]}"; do
    if [[ ! -f "$f" ]]; then
        echo "Warning: Input file not found, skipping: $f" >&2
        continue
    fi
    # Extract chromosome from filename (format: {sample}_{strand}_{chr}.bg)
    # Strand is always 'for' or 'rev', use it as delimiter since chr may contain underscores
    base=$(basename "$f")
    # Remove extension
    base_noext="${base%.*}"
    # Extract chromosome: everything after _for_ or _rev_
    if [[ "$base_noext" == *"_for_"* ]]; then
        chr="${base_noext#*_for_}"
    elif [[ "$base_noext" == *"_rev_"* ]]; then
        chr="${base_noext#*_rev_}"
    else
        echo "Warning: Cannot extract chromosome from filename (no _for_ or _rev_): $f" >&2
        continue
    fi
    CHR_NAMES+=("$chr")
    CHR_FILES+=("$f")
done

# Helper function to lookup file by chromosome
get_file_for_chr() {
    local target="$1"
    local i
    for ((i=0; i<${#CHR_NAMES[@]}; i++)); do
        if [[ "${CHR_NAMES[$i]}" == "$target" ]]; then
            echo "${CHR_FILES[$i]}"
            return 0
        fi
    done
    return 1
}

# Clear output file
> "$OUTPUT"

# Concatenate files in genome order
echo "Merging in genome order..."
merged_count=0
while read -r chr size rest; do
    file=$(get_file_for_chr "$chr") && {
        cat "$file" >> "$OUTPUT"
        echo "  + $chr ($file)"
        ((merged_count++))
    }
done < "$GENOME_FAI"

echo ""
echo "Done. Merged $merged_count files into $OUTPUT"
