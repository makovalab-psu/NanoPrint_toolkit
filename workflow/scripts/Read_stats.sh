#!/bin/bash

# Read Stats: Calculate sequencing statistics from FASTQ or BAM files
# Outputs: File prefix, Gigabases, Reads (millions), N50, Q50

#set -euo pipefail

# Usage function
usage() {
    cat << EOF
Usage: $(basename "$0") -i <input.fastq.gz|input.bam> -o <output.txt>

Calculate read statistics from FASTQ or BAM files using streaming mode.
Memory-efficient: uses O(1) memory regardless of input file size.

Required arguments:
    -i    Input file (*.fastq.gz or *.bam)
    -o    Output file (tab-delimited)

Optional arguments:
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

while getopts "i:o:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
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

# Get output directory and create it
OUT_DIR=$(dirname "$OUTPUT")
if [[ -z "$OUT_DIR" || "$OUT_DIR" == "." ]]; then
    OUT_DIR="$(pwd)"
fi
mkdir -p "$OUT_DIR"

echo "Calculating read statistics (streaming mode)..."

# Memory-efficient streaming approach:
# - Use length histogram (array indexed by length) instead of storing all values
# - Use quality histogram (binned to 0.1 precision) for Q50
# - Single pass through data, O(1) memory regardless of file size

# Extract and compute statistics in a single streaming pass
if [[ "$FILE_TYPE" == "fastq" ]]; then
    # FASTQ format: 4 lines per read (header, seq, +, qual)
    STATS=$(gunzip -c "$INPUT" | awk '
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
    }
    NR % 4 == 0 {
        # Quality line - calculate mean Phred score
        qual_str = $0
        qual_len = length(qual_str)
        sum_error = 0
        for (i = 1; i <= qual_len; i++) {
            char = substr(qual_str, i, 1)
            phred = ord[char] - 33
            sum_error += 10^(-phred/10)
        }
        mean_error = (qual_len > 0) ? sum_error / qual_len : 1
        mean_phred = -10 * log(mean_error)/log(10)

        # Accumulate into histograms
        len_hist[len] += len        # Total bases at this length
        len_count[len]++            # Count of reads at this length

        # Bin quality to 0.1 precision (multiply by 10, round)
        qual_bin = int(mean_phred * 10 + 0.5)
        qual_hist[qual_bin] += len  # Weight by read length (bases)

        total_bp += len
        read_count++
    }
    END {
        if (read_count == 0) {
            print "ERROR:0:0:0:0"
            exit
        }

        # Find N50: traverse lengths from high to low
        half_bp = total_bp / 2
        cumsum = 0
        n50 = 0
        # Find max length for iteration
        max_len = 0
        for (l in len_hist) {
            if (l + 0 > max_len) max_len = l + 0
        }
        for (l = max_len; l >= 1; l--) {
            if (l in len_hist) {
                cumsum += len_hist[l]
                if (cumsum >= half_bp) {
                    n50 = l
                    break
                }
            }
        }

        # Find Q50: traverse quality bins from high to low
        cumsum = 0
        q50 = 0
        max_qual = 0
        for (q in qual_hist) {
            if (q + 0 > max_qual) max_qual = q + 0
        }
        for (q = max_qual; q >= 0; q--) {
            if (q in qual_hist) {
                cumsum += qual_hist[q]
                if (cumsum >= half_bp) {
                    q50 = q / 10.0  # Convert back from bin
                    break
                }
            }
        }

        print total_bp ":" read_count ":" n50 ":" q50
    }
    ')
elif [[ "$FILE_TYPE" == "bam" ]]; then
    # BAM format: use samtools to extract
    STATS=$(samtools view "$INPUT" | awk '
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
        qual_len = length(qual_str)
        sum_error = 0
        for (i = 1; i <= qual_len; i++) {
            char = substr(qual_str, i, 1)
            phred = ord[char] - 33
            sum_error += 10^(-phred/10)
        }
        mean_error = (qual_len > 0) ? sum_error / qual_len : 1
        mean_phred = -10 * log(mean_error)/log(10)

        # Accumulate into histograms
        len_hist[len] += len
        len_count[len]++
        qual_bin = int(mean_phred * 10 + 0.5)
        qual_hist[qual_bin] += len

        total_bp += len
        read_count++
    }
    END {
        if (read_count == 0) {
            print "ERROR:0:0:0:0"
            exit
        }

        half_bp = total_bp / 2
        cumsum = 0
        n50 = 0
        max_len = 0
        for (l in len_hist) {
            if (l + 0 > max_len) max_len = l + 0
        }
        for (l = max_len; l >= 1; l--) {
            if (l in len_hist) {
                cumsum += len_hist[l]
                if (cumsum >= half_bp) {
                    n50 = l
                    break
                }
            }
        }

        cumsum = 0
        q50 = 0
        max_qual = 0
        for (q in qual_hist) {
            if (q + 0 > max_qual) max_qual = q + 0
        }
        for (q = max_qual; q >= 0; q--) {
            if (q in qual_hist) {
                cumsum += qual_hist[q]
                if (cumsum >= half_bp) {
                    q50 = q / 10.0
                    break
                }
            }
        }

        print total_bp ":" read_count ":" n50 ":" q50
    }
    ')
fi

# Parse results
TOTAL_BP=$(echo "$STATS" | cut -d: -f1)
READ_COUNT=$(echo "$STATS" | cut -d: -f2)
N50=$(echo "$STATS" | cut -d: -f3)
Q50=$(echo "$STATS" | cut -d: -f4)

echo "Total reads: $READ_COUNT"

if [[ "$READ_COUNT" -eq 0 || "$STATS" == ERROR* ]]; then
    echo "Error: No reads found in input file" >&2
    exit 1
fi

echo "N50: $N50"
echo "Q50: $Q50"

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


