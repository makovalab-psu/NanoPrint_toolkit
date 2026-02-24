#!/bin/bash

# Merge BigWig: Merge chromosome-split bigWig files in genome order
# Converts to bedGraph, concatenates, then converts back to bigWig

set -euo pipefail

# Usage function
usage() {
    cat << EOF
Usage: $(basename "$0") -g <genome.fa.fai> -o <output.bw> [-T tmpdir] <input1.bw> <input2.bw> ...

Merge chromosome-split bigWig files in genome order.

Required arguments:
    -g    Genome sizes file (samtools faidx .fai format)
    -o    Output bigWig file

Optional arguments:
    -T    Temporary directory (default: same directory as output)
    -h    Show this help message

Positional arguments:
    Remaining arguments are input bigWig files to merge

The input files are merged in the order of chromosomes as they appear
in the genome sizes file (first column of .fai).

Dependencies:
    - bigWigToBedGraph (UCSC tools)
    - bedGraphToBigWig (UCSC tools)

Example:
    $(basename "$0") -g genome.fa.fai -o merged.bw sample_chr1.bw sample_chr2.bw sample_chrX.bw
EOF
    exit 1
}

# Parse arguments
GENOME_FAI=""
OUTPUT=""
TMP_DIR=""

while getopts "g:o:T:h" opt; do
    case $opt in
        g) GENOME_FAI="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        T) TMP_DIR="$OPTARG" ;;
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

# Check required tools
if ! command -v bigWigToBedGraph &> /dev/null; then
    echo "Error: bigWigToBedGraph not found. Please install UCSC tools." >&2
    exit 1
fi

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
    TMP_DIR="${OUT_DIR}/tmp_merge_bw_$$"
fi
mkdir -p "$TMP_DIR"

# Cleanup function
cleanup() {
    if [[ -d "$TMP_DIR" ]]; then
        rm -rf "$TMP_DIR"
    fi
}
trap cleanup EXIT INT TERM

echo "=== Merge BigWig Files ==="
echo "Genome: $GENOME_FAI"
echo "Output: $OUTPUT"
echo "Input files: ${#INPUT_FILES[@]}"
echo "Temp dir: $TMP_DIR"
echo ""

# Create chrom.sizes file from fai (columns 1 and 2)
CHROM_SIZES="${TMP_DIR}/chrom.sizes"
cut -f1,2 "$GENOME_FAI" > "$CHROM_SIZES"

# Parallel arrays for Bash 3.2 compatibility (macOS)
CHR_NAMES=()
CHR_FILES=()

# Lookup function to get file for a chromosome
get_file_for_chr() {
    local target_chr="$1"
    local i
    for ((i=0; i<${#CHR_NAMES[@]}; i++)); do
        if [[ "${CHR_NAMES[i]}" == "$target_chr" ]]; then
            echo "${CHR_FILES[i]}"
            return 0
        fi
    done
    return 1
}

# Build arrays of input files by chromosome
for f in "${INPUT_FILES[@]}"; do
    if [[ ! -f "$f" ]]; then
        echo "Warning: Input file not found, skipping: $f" >&2
        continue
    fi
    # Extract chromosome from filename (format: {sample}_{strand}_{chr}.bw)
    # Strand is always 'for' or 'rev', use it as delimiter since chr may contain underscores
    base=$(basename "$f")
    base_noext="${base%.bw}"
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

# Convert each bigWig to bedGraph and concatenate in genome order
echo "Converting bigWig files to bedGraph and merging..."
MERGED_BG="${TMP_DIR}/merged.bg"
> "$MERGED_BG"

merged_count=0
while read -r chr size rest; do
    bw_file=$(get_file_for_chr "$chr")
    if [[ -n "$bw_file" ]]; then
        tmp_bg="${TMP_DIR}/${chr}.bg"

        echo "  Converting $chr..."
        bigWigToBedGraph "$bw_file" "$tmp_bg"
        cat "$tmp_bg" >> "$MERGED_BG"
        rm "$tmp_bg"

        merged_count=$((merged_count + 1))
    fi
done < "$GENOME_FAI"

echo ""
echo "Merged $merged_count chromosomes"

# Sort bedGraph (required by bedGraphToBigWig)
echo "Sorting merged bedGraph..."
SORTED_BG="${TMP_DIR}/merged_sorted.bg"
sort -k1,1 -k2,2n "$MERGED_BG" > "$SORTED_BG"

# Convert to bigWig
echo "Converting to bigWig..."
bedGraphToBigWig "$SORTED_BG" "$CHROM_SIZES" "$OUTPUT"

echo ""
echo "Done. Output: $OUTPUT"
