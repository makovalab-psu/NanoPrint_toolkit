#!/bin/bash

# Read Stats: Calculate sequencing statistics from FASTQ or BAM files
# Outputs: File prefix, Gigabases, Reads (millions), N50, Q50

#set -euo pipefail

# Usage function
usage() {
    cat << EOF
Usage: $(basename "$0") -i <input.fastq.gz|input.bam> -o <output.txt> [-T tmpdir]

Calculate read statistics from FASTQ or BAM files.

Required arguments:
    -i    Input file (*.fastq.gz or *.bam)
    -o    Output file (tab-delimited)

Optional arguments:
    -T    Temporary directory (default: same directory as output)
    -h    Show this help message

Input formats:
    - Gzipped FASTQ (*.fastq.gz)
    - Unmapped BAM (*.bam)

Output format (tab-delimited with header):
    1. File_prefix    - Input filename without extension
    2. Giga_bp        - Total gigabases
    3. Reads_million  - Total reads in millions
    4. N50            - Read length N50 (length at 50% of total bases)
    5. Q50            - Read quality Q50 (mean Phred at 50% of total bases)

Example:
    $(basename "$0") -i sample.fastq.gz -o sample_stats.txt
    $(basename "$0") -i sample.bam -o sample_stats.txt
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

# Determine file type from extension
FILENAME=$(basename "$INPUT")
if [[ "$FILENAME" == *.fastq.gz ]]; then
    FILE_TYPE="fastq"
    FILE_PREFIX="${FILENAME%.fastq.gz}"
elif [[ "$FILENAME" == *.fq.gz ]]; then
    FILE_TYPE="fastq"
    FILE_PREFIX="${FILENAME%.fq.gz}"
elif [[ "$FILENAME" == *.bam ]]; then
    FILE_TYPE="bam"
    FILE_PREFIX="${FILENAME%.bam}"
    # Check samtools is available
    if ! command -v samtools &> /dev/null; then
        echo "Error: samtools not found. Required for BAM input." >&2
        exit 1
    fi
else
    echo "Error: Unrecognized file type. Expected *.fastq.gz or *.bam" >&2
    exit 1
fi

echo "Input file: $INPUT"
echo "File type: $FILE_TYPE"
echo "File prefix: $FILE_PREFIX"

# Get output directory
OUT_DIR=$(dirname "$OUTPUT")
if [[ -z "$OUT_DIR" || "$OUT_DIR" == "." ]]; then
    OUT_DIR="$(pwd)"
fi
mkdir -p "$OUT_DIR"

# Set temp directory (default: in output directory)
if [[ -z "$TMP_DIR" ]]; then
    TMP_DIR="${OUT_DIR}/tmp_stats_$$"
fi
mkdir -p "$TMP_DIR"

# Temp file for length and quality data
LENGTH_QUAL="${TMP_DIR}/length_qual.txt"

echo "Extracting read lengths and qualities..."

# Extract sequence lengths and mean quality scores
if [[ "$FILE_TYPE" == "fastq" ]]; then
    # FASTQ format: 4 lines per read (header, seq, +, qual)
    gunzip -c "$INPUT" | awk '
    BEGIN {
        OFS = "\t"
        # Create lookup array for ASCII to integer conversion
        for (n = 0; n < 256; n++) {
            ord[sprintf("%c", n)] = n
        }
    }
    NR % 4 == 2 { 
        # Sequence line - get length
        len = length($0)
        seq = $0
    }
    NR % 4 == 0 {
        # Quality line - calculate mean Phred score
        qual_str = $0
        sum_phred = 0
        qual_len = length(qual_str)
        for (i = 1; i <= qual_len; i++) {
            char = substr(qual_str, i, 1)
            phred = ord[char] - 33
            sum_error += 10^(-phred/10)
        }
        mean_error = (qual_len > 0) ?  sum_error / qual_len : 0
        mean_phred = -10 * log(mean_error)/log(10)
        print len, mean_phred
    }
    ' > "$LENGTH_QUAL"
elif [[ "$FILE_TYPE" == "bam" ]]; then
    # BAM format: use samtools to extract
    samtools view "$INPUT" | awk '
    BEGIN {
        OFS = "\t"
        # Create lookup array for ASCII to integer conversion
        for (n = 0; n < 256; n++) {
            ord[sprintf("%c", n)] = n
        }
    }
    {
        # Column 10 = sequence, Column 11 = quality
        len = length($10)
        qual_str = $11
        sum_phred = 0
        qual_len = length(qual_str)
        for (i = 1; i <= qual_len; i++) {
            char = substr(qual_str, i, 1)
            phred = ord[char] - 33
            sum_error += 10^(-phred/10)
        }
        mean_error = (qual_len > 0) ?  sum_error / qual_len : 0
        mean_phred = -10 * log(mean_error)/log(10)
        print len, mean_phred
    }
    ' > "$LENGTH_QUAL"
fi

READ_COUNT=$(wc -l < "$LENGTH_QUAL" | tr -d ' ')
echo "Total reads: $READ_COUNT"

if [[ "$READ_COUNT" -eq 0 ]]; then
    echo "Error: No reads found in input file" >&2
    exit 1
fi

echo "Calculating N50..."
# Calculate total bases and N50
# Sort by length descending, accumulate until 50% of total bases
TOTAL_BP=$(awk '{sum += $1} END {print sum}' "$LENGTH_QUAL")
N50=$(sort -t$'\t' -k1,1 -rn "$LENGTH_QUAL" | awk -v total="$TOTAL_BP" '
BEGIN { cumsum = 0; half = total / 2 }
{
    cumsum += $1
    if (cumsum >= half) {
        print $1
        exit
    }
}
')

echo "Calculating Q50..."
# Calculate Q50: sort by mean quality descending, find quality at 50% of bases
Q50=$(sort -t$'\t' -k2,2 -rn "$LENGTH_QUAL" | awk -v total="$LENGTH_QUAL" '
BEGIN { cumsum = 0; half = total / 2 }
{
    cumsum += 1
    if (cumsum >= half) {
        print $2
        exit
    }
}
')

# Calculate summary statistics
GIGA_BP=$(echo "scale=4; $TOTAL_BP / 1000000000" | bc)
READS_MILLION=$(echo "scale=4; $READ_COUNT / 1000000" | bc)

# Write output
{
    echo -e "File_prefix\tGiga_bp\tReads_million\tN50\tQ50"
    printf "%s\t%s\t%s\t%s\t%.2f\n" "$FILE_PREFIX" "$GIGA_BP" "$READS_MILLION" "$N50" "$Q50"
} > "$OUTPUT"

# Display results
echo ""
echo "=== Read Statistics ==="
cat "$OUTPUT" | column -t
echo ""
echo "Done: $OUTPUT"

# Cleanup function
cleanup() {
    if [[ -d "$TMP_DIR" ]]; then
        rm -rf "$TMP_DIR"
    fi
}
trap cleanup EXIT


