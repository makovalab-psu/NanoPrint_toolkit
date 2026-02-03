#!/bin/bash

# Calculate Reactivity: Treatment (MnO4) minus Control
# Simplified script for permanganate footprinting analysis

set -euo pipefail

# Get script directory for finding Python scripts
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default coverage threshold
COV_THRESHOLD=10

# Usage function
usage() {
    cat << EOF
Usage: $(basename "$0") -p <MnO4.txt.gz> -m <CTRL.txt.gz> -o <output.txt.gz> [-c threshold] [-T tmpdir]

Calculate reactivity from perbase error (treatment minus control).

Required arguments:
    -p    Treatment file (MnO4/plus, gzipped)
    -m    Control file (CTRL/minus, gzipped)
    -o    Output file (gzipped)

Optional arguments:
    -c    Minimum coverage threshold (default: 10)
    -T    Temporary directory (default: same directory as output)
    -h    Show this help message

Input format (5 columns, tab-separated):
    1. Chromosome name
    2. Position (1-based)
    3. Nucleotide identity
    4. Coverage
    5. Perbase error

Output format (4 columns, tab-separated):
    1. Chromosome name
    2. Position (1-based)
    3. Nucleotide identity
    4. Reactivity (treatment_error - control_error)
       Special values:
         999999  = position missing in control
        -999999  = position missing in treatment

Notes:
    - Input files must contain exactly one chromosome each
    - Chromosome names must match between treatment and control files
    - Chromosome size is determined from the maximum position in input files

Example:
    $(basename "$0") -p MnO4_chr1.txt.gz -m CTRL_chr1.txt.gz -o react_chr1.txt.gz
    $(basename "$0") -p MnO4_chr1.txt.gz -m CTRL_chr1.txt.gz -o react_chr1.txt.gz -c 20 -T /tmp
EOF
    exit 1
}

# Parse arguments
PLUS=""
MINUS=""
OUT=""
TMP_DIR=""

while getopts "p:m:o:c:T:h" opt; do
    case $opt in
        p) PLUS="$OPTARG" ;;
        m) MINUS="$OPTARG" ;;
        o) OUT="$OPTARG" ;;
        c) COV_THRESHOLD="$OPTARG" ;;
        T) TMP_DIR="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Validate required arguments
if [[ -z "$PLUS" || -z "$MINUS" || -z "$OUT" ]]; then
    echo "Error: Missing required arguments" >&2
    usage
fi

# Check input files exist
if [[ ! -f "$PLUS" ]]; then
    echo "Error: Treatment file not found: $PLUS" >&2
    exit 1
fi

if [[ ! -f "$MINUS" ]]; then
    echo "Error: Control file not found: $MINUS" >&2
    exit 1
fi

# Check Python script exists
FILL_SCRIPT="${SCRIPT_DIR}/Fill_in_blanks_for_perbase_error_files.py"
if [[ ! -f "$FILL_SCRIPT" ]]; then
    echo "Error: Fill_in_blanks script not found: $FILL_SCRIPT" >&2
    exit 1
fi

# Create output directory if needed
OUT_DIR=$(dirname "$OUT")
if [[ -n "$OUT_DIR" && "$OUT_DIR" != "." ]]; then
    mkdir -p "$OUT_DIR"
fi

# Set temp directory (default: same as output directory)
if [[ -z "$TMP_DIR" ]]; then
    TMP_DIR="${OUT_DIR:-.}/tmp_react_$$"
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
PLUS_DECOMP="${TMP_DIR}/plus_decomp.txt"
MINUS_DECOMP="${TMP_DIR}/minus_decomp.txt"
PLUS_FILLED="${TMP_DIR}/plus_filled.txt"
MINUS_FILLED="${TMP_DIR}/minus_filled.txt"

# Decompress input files (using gunzip -c for macOS/Linux compatibility)
echo "Decompressing treatment file..."
gunzip -c "$PLUS" > "$PLUS_DECOMP"

echo "Decompressing control file..."
gunzip -c "$MINUS" > "$MINUS_DECOMP"

# Get chromosome name from first line of treatment file
CHR=$(head -n 1 "$PLUS_DECOMP" | cut -f1)
echo "Chromosome: $CHR"

# Determine chromosome size from max position in both files
echo "Determining chromosome size..."
PLUS_MAX=$(tail -n 1 "$PLUS_DECOMP" | cut -f2)
MINUS_MAX=$(tail -n 1 "$MINUS_DECOMP" | cut -f2)

if [[ "$PLUS_MAX" -gt "$MINUS_MAX" ]]; then
    CHR_SIZE="$PLUS_MAX"
else
    CHR_SIZE="$MINUS_MAX"
fi

echo "Chromosome size: $CHR_SIZE"

# Fill in missing positions to align both files
echo "Filling missing positions in treatment file..."
python "$FILL_SCRIPT" "$PLUS_DECOMP" "$CHR_SIZE" "$PLUS_FILLED"

echo "Filling missing positions in control file..."
python "$FILL_SCRIPT" "$MINUS_DECOMP" "$CHR_SIZE" "$MINUS_FILLED"

# Process files using awk
# Join on chr+position, calculate reactivity
echo "Calculating reactivity..."
paste "$PLUS_FILLED" "$MINUS_FILLED" | awk -v cov="$COV_THRESHOLD" '
BEGIN { FS="\t"; OFS="\t" }
{
    # Treatment (MnO4): columns 1-5
    # Control (CTRL): columns 6-10
    chr_t = $1; pos_t = $2; nuc_t = $3; cov_t = $4; err_t = $5
    chr_c = $6; pos_c = $7; nuc_c = $8; cov_c = $9; err_c = $10

    # Check if position exists in both files (non-empty, has all columns)
    has_treatment = (chr_t != "" && nuc_t != "" && err_t != "")
    has_control = (chr_c != "" && nuc_c != "" && err_c != "")

    # Determine chromosome, position, nucleotide from available data
    if (has_treatment) {
        chr = chr_t; pos = pos_t; nuc = nuc_t
    } else if (has_control) {
        chr = chr_c; pos = pos_c; nuc = nuc_c
    } else {
        next  # Skip if both are empty/placeholder lines
    }

    # Calculate reactivity based on data availability
    if (has_treatment && has_control) {
        # Both present - check coverage threshold
        if (cov_t >= cov && cov_c >= cov) {
            reactivity = err_t - err_c
            print chr, pos, nuc, reactivity
        }
    } else if (has_treatment && !has_control) {
        # Missing in control
        print chr, pos, nuc, 999999
    } else if (!has_treatment && has_control) {
        # Missing in treatment
        print chr, pos, nuc, -999999
    }
}
' | gzip > "$OUT"

echo "Done: $OUT"
