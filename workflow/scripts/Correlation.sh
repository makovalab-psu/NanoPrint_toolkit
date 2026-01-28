#!/bin/bash

# Correlation: Correlate per-base error between two samples via subsampling
# Randomly samples genomic positions and compares error rates between samples

set -euo pipefail

# Usage function
usage() {
    cat << EOF
Usage: $(basename "$0") -a <perbase_error1.txt.gz> -b <perbase_error2.txt.gz> -g <genome.fa.fai> -o <output.txt> [-s subsamples] [-T tmpdir]

Correlate per-base error rates between two samples by random subsampling.

Required arguments:
    -a    Per-base error file 1 (gzipped)
    -b    Per-base error file 2 (gzipped)
    -g    Genome sizes file (samtools faidx .fai format)
    -o    Output correlation file (tab-delimited)

Optional arguments:
    -s    Number of subsamples (default: 10000)
    -T    Temporary directory (default: same directory as output)
    -h    Show this help message

Input format (per-base error files, gzipped):
    Column 1: Chromosome
    Column 2: Position (1-based)
    Column 3: Nucleotide
    Column 4: Coverage
    Column 5: Per-base error

Output format (tab-delimited with header):
    Column 1: Chromosome
    Column 2: Nucleotide
    Column 3: Perbase_error_1
    Column 4: Perbase_error_2
    Column 5: Coverage_1
    Column 6: Coverage_2

Example:
    $(basename "$0") -a sample1.txt.gz -b sample2.txt.gz -g genome.fa.fai -o correlation.txt -s 10000
EOF
    exit 1
}

# Parse arguments
FILE_A=""
FILE_B=""
GENOME_FAI=""
OUTPUT=""
SUBSAMPLES=10000
TMP_DIR=""

while getopts "a:b:g:o:s:T:h" opt; do
    case $opt in
        a) FILE_A="$OPTARG" ;;
        b) FILE_B="$OPTARG" ;;
        g) GENOME_FAI="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        s) SUBSAMPLES="$OPTARG" ;;
        T) TMP_DIR="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Validate required arguments
if [[ -z "$FILE_A" || -z "$FILE_B" || -z "$GENOME_FAI" || -z "$OUTPUT" ]]; then
    echo "Error: Missing required arguments" >&2
    usage
fi

# Check input files exist
if [[ ! -f "$FILE_A" ]]; then
    echo "Error: Per-base error file 1 not found: $FILE_A" >&2
    exit 1
fi

if [[ ! -f "$FILE_B" ]]; then
    echo "Error: Per-base error file 2 not found: $FILE_B" >&2
    exit 1
fi

if [[ ! -f "$GENOME_FAI" ]]; then
    echo "Error: Genome sizes file not found: $GENOME_FAI" >&2
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
    TMP_DIR="${OUT_DIR}/tmp_correlation_$$"
fi
mkdir -p "$TMP_DIR"

# Cleanup function
cleanup() {
    if [[ -d "$TMP_DIR" ]]; then
        rm -rf "$TMP_DIR"
    fi
}
trap cleanup EXIT INT TERM

echo "=== Per-base Error Correlation ==="
echo "File A: $FILE_A"
echo "File B: $FILE_B"
echo "Genome: $GENOME_FAI"
echo "Subsamples: $SUBSAMPLES"
echo "Output: $OUTPUT"
echo ""

# Temp files
RANDOM_POS="${TMP_DIR}/random_positions.txt"
SORTED_POS="${TMP_DIR}/sorted_positions.txt"
DATA_A="${TMP_DIR}/data_a.txt"
DATA_B="${TMP_DIR}/data_b.txt"
MERGED="${TMP_DIR}/merged.txt"

# Step 1: Calculate total genome size and generate random positions
echo "Generating $SUBSAMPLES random genomic positions..."

# Read chromosome sizes and calculate cumulative positions for random sampling
awk -v n="$SUBSAMPLES" -v seed="$RANDOM" '
BEGIN {
    srand(seed)
    total = 0
}
{
    chr[NR] = $1
    size[NR] = $2
    cumsum[NR] = total + $2
    total += $2
    nchr = NR
}
END {
    # Generate n random positions
    for (i = 1; i <= n; i++) {
        # Random position in genome
        pos = int(rand() * total) + 1

        # Find which chromosome this falls in
        for (j = 1; j <= nchr; j++) {
            if (pos <= cumsum[j]) {
                if (j == 1) {
                    chr_pos = pos
                } else {
                    chr_pos = pos - cumsum[j-1]
                }
                print chr[j] "\t" chr_pos
                break
            }
        }
    }
}
' "$GENOME_FAI" > "$RANDOM_POS"

# Sort positions for efficient lookup
sort -k1,1 -k2,2n "$RANDOM_POS" > "$SORTED_POS"

echo "  Generated random positions"

# Step 2: Extract data from file A at sampled positions
echo "Extracting data from file A..."
gunzip -c "$FILE_A" | sort -k1,1 -k2,2n | join -t$'\t' -1 1 -2 1 \
    -o '1.1,1.2,2.2,1.3,1.4,1.5' \
    - <(awk '{print $1":"$2"\t"$2}' "$SORTED_POS" | sort -k1,1) \
    2>/dev/null | awk -F'\t' '{split($1,a,":"); print a[1]"\t"$3"\t"$4"\t"$5"\t"$6}' > "$DATA_A" || true

# Alternative approach: use awk to match positions
gunzip -c "$FILE_A" | awk -F'\t' '
    NR==FNR {positions[$1":"$2]=1; next}
    ($1":"$2) in positions {print $1"\t"$2"\t"$3"\t"$4"\t"$5}
' "$SORTED_POS" - > "$DATA_A"


# Step 3: Extract data from file B at sampled positions
echo "Extracting data from file B..."
gunzip -c "$FILE_B" | awk -F'\t' '
    NR==FNR {positions[$1":"$2]=1; next}
    ($1":"$2) in positions {print $1"\t"$2"\t"$3"\t"$4"\t"$5}
' "$SORTED_POS" - > "$DATA_B"


# Step 4: Merge data from both files
echo "Merging data..."

# Create lookup from file A and B, then merge on sampled positions
awk -F'\t' '
BEGIN {
    OFS = "\t"
    print "Chromosome\tNucleotide\tPerbase_error_1\tPerbase_error_2\tCoverage_1\tCoverage_2"
}
# Read file A data
ARGIND==1 {
    key = $1":"$2
    nuc_a[key] = $3
    cov_a[key] = $4
    err_a[key] = $5
    next
}
# Read file B data
ARGIND==2 {
    key = $1":"$2
    nuc_b[key] = $3
    cov_b[key] = $4
    err_b[key] = $5
    next
}
# Read sampled positions and output merged data
ARGIND==3 {
    key = $1":"$2
    chr = $1

    # Get nucleotide (prefer A, fall back to B, or empty)
    if (key in nuc_a) {
        nuc = nuc_a[key]
    } else if (key in nuc_b) {
        nuc = nuc_b[key]
    } else {
        nuc = ""
    }

    # Get values or empty string if missing
    e1 = (key in err_a) ? err_a[key] : ""
    e2 = (key in err_b) ? err_b[key] : ""
    c1 = (key in cov_a) ? cov_a[key] : ""
    c2 = (key in cov_b) ? cov_b[key] : ""

    print chr, nuc, e1, e2, c1, c2
}
' "$DATA_A" "$DATA_B" "$SORTED_POS" > "$OUTPUT"

echo ""
echo "Done. Output: $OUTPUT"
