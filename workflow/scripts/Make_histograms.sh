#!/bin/bash

# Make Histograms: Extract histogram data from samtools stats output
# Extracts read length, mapping quality, insertions, deletions, coverage, and GC coverage

set -euo pipefail

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Usage function
usage() {
    cat << EOF
Usage: $(basename "$0") -i <input_stats.txt> -o <output_histograms.txt> [-p <output_histograms.pdf>]

Extract histogram data from samtools stats output files.

Required arguments:
    -i    Input stats file (from samtools stats / Alignment_stats.sh)
    -o    Output histogram file (tab-delimited)

Optional arguments:
    -p    Output PDF file with histogram plots
    -h    Show this help message

Output format (tab-delimited with header):
    Column 1: Var    - Histogram type (RL, MAPQ, INS, DEL, COV)
    Column 2: Value  - Bin value (length, quality, size, coverage, GC%)
    Column 3: Count  - Count for that bin

Histogram types extracted:
    RL   - Read length distribution
    MAPQ - Mapping quality distribution
    INS  - Insertion size distribution
    DEL  - Deletion size distribution
    COV  - Coverage distribution

Example:
    $(basename "$0") -i sample_stats.txt -o sample_histograms.txt
    $(basename "$0") -i sample_stats.txt -o sample_histograms.txt -p sample_histograms.pdf
EOF
    exit 1
}

# Parse arguments
INPUT=""
OUTPUT=""
PDF_OUTPUT=""

while getopts "i:o:p:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        p) PDF_OUTPUT="$OPTARG" ;;
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

# Check R script exists if PDF output requested
R_SCRIPT="${SCRIPT_DIR}/Plot_histograms.R"
if [[ -n "$PDF_OUTPUT" ]]; then
    if [[ ! -f "$R_SCRIPT" ]]; then
        echo "Error: R plotting script not found: $R_SCRIPT" >&2
        exit 1
    fi
    if ! command -v Rscript &> /dev/null; then
        echo "Error: Rscript not found. Please install R." >&2
        exit 1
    fi
fi

# Create output directory if needed
OUT_DIR=$(dirname "$OUTPUT")
if [[ -n "$OUT_DIR" && "$OUT_DIR" != "." ]]; then
    mkdir -p "$OUT_DIR"
fi

# Create PDF output directory if needed
if [[ -n "$PDF_OUTPUT" ]]; then
    PDF_DIR=$(dirname "$PDF_OUTPUT")
    if [[ -n "$PDF_DIR" && "$PDF_DIR" != "." ]]; then
        mkdir -p "$PDF_DIR"
    fi
fi

echo "=== Histogram Extraction ==="
echo "Input: $INPUT"
echo "Output: $OUTPUT"
if [[ -n "$PDF_OUTPUT" ]]; then
    echo "PDF Output: $PDF_OUTPUT"
fi
echo ""

# Write header
echo -e "Var\tValue\tCount" > "$OUTPUT"

# Extract Read Length histogram (RL)
# Format: RL	read_length	count
echo "Extracting read length distribution (RL)..."
grep "^RL" "$INPUT" | awk -F'\t' '{print "RL\t" $2 "\t" $3}' >> "$OUTPUT"

# Extract Mapping Quality histogram (MAPQ)
# Format: MAPQ	mapq_value	count
echo "Extracting mapping quality distribution (MAPQ)..."
grep "^MAPQ" "$INPUT" | awk -F'\t' '{print "MAPQ\t" $2 "\t" $3}' >> "$OUTPUT"

# Extract Insertion size histogram (from ID lines)
# Format: ID	length	insertions	deletions
echo "Extracting insertion size distribution (INS)..."
grep "^ID" "$INPUT" | awk -F'\t' '{print "INS\t" $2 "\t" $3}' >> "$OUTPUT"

# Extract Deletion size histogram (from ID lines)
# Format: ID	length	insertions	deletions
echo "Extracting deletion size distribution (DEL)..."
grep "^ID" "$INPUT" | awk -F'\t' '{print "DEL\t" $2 "\t" $4}' >> "$OUTPUT"

# Extract Coverage histogram (COV)
# Format: COV	[start-end]	count
echo "Extracting coverage distribution (COV)..."
grep "^COV" "$INPUT" | awk -F'\t' '{print "COV\t" $2 "\t" $3}' >> "$OUTPUT"

echo ""
echo "Histogram data: $OUTPUT"

# Generate PDF if requested
if [[ -n "$PDF_OUTPUT" ]]; then
    echo ""
    echo "Generating PDF plots..."
    Rscript "$R_SCRIPT" "$OUTPUT" "$PDF_OUTPUT"
    echo "PDF output: $PDF_OUTPUT"
fi

echo ""
echo "Done."
