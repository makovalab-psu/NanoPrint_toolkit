#!/bin/bash

# Alignment Stats: Generate summary statistics from raw and filtered BAM files
# Uses samtools stats and flagstat to extract alignment metrics

set -euo pipefail

# Usage function
usage() {
    cat << EOF
Usage: $(basename "$0") [-a <raw.bam>] -f <filtered.bam> -o <output.txt> [-d <stats_dir>]

Generate alignment statistics from raw and filtered BAM files using samtools.

Required arguments:
    -f    Filtered alignment BAM file
    -o    Output table file (tab-delimited)

Optional arguments:
    -a    Raw (pre-filter) alignment BAM file. Omit it when no pre-filter BAM
          exists — an input supplied already Uncalled4-aligned (^u in CONFIG)
          has one alignment, not two. The Raw_alignment column then reads NA.
    -d    Directory for the per-BAM stats/flagstat files (default: next to each
          BAM). Needed when the BAM lives outside the data/ tree it is reported
          under, so the files land where the histograms rule looks for them.
    -h    Show this help message

Outputs:
    For each BAM file given:
        <bam>.bai           - BAM index file (written next to the BAM)
        <name>_stats.txt     - samtools stats output
        <name>_flagstats.txt - samtools flagstat output

    Output table (tab-delimited with header):
        Column 1: Sample              - Sample name from BAM filename
        Column 2: Statistic           - Description of the statistic
        Column 3: Raw_alignment       - Value from raw alignment, or NA
        Column 4: Filtered_alignment  - Value from filtered alignment

Example:
    $(basename "$0") -a sample_raw.bam -f sample_filtered.bam -o sample_alignment_stats.txt
    $(basename "$0") -f uncalled4.bam -d data/filtered_alignments/genome -o stats.txt
EOF
    exit 1
}

# Parse arguments
RAW_BAM=""
FILTERED_BAM=""
OUTPUT=""
STATS_DIR=""

while getopts "a:f:o:d:h" opt; do
    case $opt in
        a) RAW_BAM="$OPTARG" ;;
        f) FILTERED_BAM="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        d) STATS_DIR="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Validate required arguments
if [[ -z "$FILTERED_BAM" || -z "$OUTPUT" ]]; then
    echo "Error: Missing required arguments" >&2
    usage
fi

# A raw BAM is optional; without one the raw column is reported as NA.
HAVE_RAW="true"
[[ -z "$RAW_BAM" ]] && HAVE_RAW="false"

# Check input files exist
if [[ "$HAVE_RAW" == "true" && ! -f "$RAW_BAM" ]]; then
    echo "Error: Raw BAM file not found: $RAW_BAM" >&2
    exit 1
fi

if [[ ! -f "$FILTERED_BAM" ]]; then
    echo "Error: Filtered BAM file not found: $FILTERED_BAM" >&2
    exit 1
fi

# Check samtools is available
if ! command -v samtools &> /dev/null; then
    echo "Error: samtools not found. Please install samtools." >&2
    exit 1
fi

# Sample name comes from the raw BAM when there is one, so existing outputs are
# unchanged; otherwise from the filtered BAM.
if [[ "$HAVE_RAW" == "true" ]]; then
    SAMPLE_NAME=$(basename "$RAW_BAM" .bam)
else
    SAMPLE_NAME=$(basename "$FILTERED_BAM" .bam)
fi

# Define output file paths. -d overrides the default of writing next to each BAM.
RAW_DIR="${STATS_DIR:-$(dirname "$RAW_BAM")}"
FILTERED_DIR="${STATS_DIR:-$(dirname "$FILTERED_BAM")}"
[[ -n "$STATS_DIR" ]] && mkdir -p "$STATS_DIR"

RAW_INDEX="${RAW_BAM}.bai"
RAW_STATS="${RAW_DIR}/$(basename "$RAW_BAM" .bam)_stats.txt"
RAW_FLAGSTATS="${RAW_DIR}/$(basename "$RAW_BAM" .bam)_flagstats.txt"

FILTERED_INDEX="${FILTERED_BAM}.bai"
FILTERED_STATS="${FILTERED_DIR}/$(basename "$FILTERED_BAM" .bam)_stats.txt"
FILTERED_FLAGSTATS="${FILTERED_DIR}/$(basename "$FILTERED_BAM" .bam)_flagstats.txt"

echo "=== Alignment Statistics Pipeline ==="
echo "Raw BAM: ${RAW_BAM:-none (filtered only)}"
echo "Filtered BAM: $FILTERED_BAM"
echo "Sample name: $SAMPLE_NAME"
echo ""

# Step 1: Index BAM files
echo "Indexing BAM files..."
if [[ "$HAVE_RAW" == "true" ]]; then
    samtools index "$RAW_BAM"
    echo "  Created: $RAW_INDEX"
fi
[[ -f "$FILTERED_INDEX" ]] || samtools index "$FILTERED_BAM"
echo "  Created: $FILTERED_INDEX"

# Step 2: Generate stats files
echo "Generating samtools stats..."
if [[ "$HAVE_RAW" == "true" ]]; then
    samtools stats "$RAW_BAM" > "$RAW_STATS"
    echo "  Created: $RAW_STATS"
fi
samtools stats "$FILTERED_BAM" > "$FILTERED_STATS"
echo "  Created: $FILTERED_STATS"

# Step 3: Generate flagstats files
echo "Generating samtools flagstat..."
if [[ "$HAVE_RAW" == "true" ]]; then
    samtools flagstat "$RAW_BAM" > "$RAW_FLAGSTATS"
    echo "  Created: $RAW_FLAGSTATS"
fi
samtools flagstat "$FILTERED_BAM" > "$FILTERED_FLAGSTATS"
echo "  Created: $FILTERED_FLAGSTATS"

# Step 4: Parse stats and flagstats into tidy output table
echo "Parsing statistics into output table..."

# Create output directory if needed
OUT_DIR=$(dirname "$OUTPUT")
if [[ -n "$OUT_DIR" && "$OUT_DIR" != "." ]]; then
    mkdir -p "$OUT_DIR"
fi

# Function to extract value from samtools stats (SN lines)
# Format: "SN	statistic_name:	value	# comment"
extract_stat() {
    local stats_file="$1"
    local pattern="$2"
    grep "^SN" "$stats_file" | grep "$pattern" | awk -F'\t' '{print $3}'
}

# Function to extract value from flagstats
# Format: "12345 + 0 description"
extract_flagstat() {
    local flagstats_file="$1"
    local pattern="$2"
    grep "$pattern" "$flagstats_file" | awk '{print $1}'
}

# Write header
echo -e "Sample\tStatistic\tRaw_alignment\tFiltered_alignment" > "$OUTPUT"

# Extract and write stats from samtools stats
# Using associative array for stat descriptions and grep patterns
declare -a STATS_PATTERNS=(
    "raw total sequences:"
    "total length:"
    "bases mapped:"
    "bases mapped (cigar):"
    "mismatches:"
    "error rate:"
    "average length:"
    "average quality:"
)

declare -a STATS_NAMES=(
    "Total sequences"
    "Total length"
    "Bases mapped"
    "Bases mapped (cigar)"
    "Mismatches"
    "Error rate"
    "Average length"
    "Average quality"
)

for i in "${!STATS_PATTERNS[@]}"; do
    pattern="${STATS_PATTERNS[$i]}"
    name="${STATS_NAMES[$i]}"
    if [[ "$HAVE_RAW" == "true" ]]; then
        raw_val=$(extract_stat "$RAW_STATS" "$pattern")
    else
        raw_val="NA"
    fi
    filtered_val=$(extract_stat "$FILTERED_STATS" "$pattern")
    echo -e "${SAMPLE_NAME}\t${name}\t${raw_val}\t${filtered_val}" >> "$OUTPUT"
done

# Extract and write stats from flagstats
# Flagstat patterns and their descriptions
declare -a FLAGSTAT_PATTERNS=(
    "in total"
    "primary$"
    "secondary$"
    "supplementary$"
    "duplicates$"
    "primary duplicates"
    "^[0-9].*mapped \("
    "primary mapped"
    "paired in sequencing"
    "read1$"
    "read2$"
    "properly paired"
    "with itself and mate mapped"
    "singletons"
    "with mate mapped to a different chr$"
    "with mate mapped to a different chr \(mapQ>=5\)"
)

declare -a FLAGSTAT_NAMES=(
    "Total reads (QC-passed + QC-failed)"
    "Primary alignments"
    "Secondary alignments"
    "Supplementary alignments"
    "Duplicates"
    "Primary duplicates"
    "Mapped reads"
    "Primary mapped"
    "Paired in sequencing"
    "Read1"
    "Read2"
    "Properly paired"
    "With itself and mate mapped"
    "Singletons"
    "Mate mapped to different chr"
    "Mate mapped to different chr (mapQ>=5)"
)

for i in "${!FLAGSTAT_PATTERNS[@]}"; do
    pattern="${FLAGSTAT_PATTERNS[$i]}"
    name="${FLAGSTAT_NAMES[$i]}"
    if [[ "$HAVE_RAW" == "true" ]]; then
        raw_val=$(grep -E "$pattern" "$RAW_FLAGSTATS" | head -1 | awk '{print $1}')
    else
        raw_val="NA"
    fi
    filtered_val=$(grep -E "$pattern" "$FILTERED_FLAGSTATS" | head -1 | awk '{print $1}')
    # Only write if we got values (pattern matched)
    if [[ -n "$raw_val" && -n "$filtered_val" ]]; then
        echo -e "${SAMPLE_NAME}\t${name}\t${raw_val}\t${filtered_val}" >> "$OUTPUT"
    fi
done

# Display results
echo ""
echo "=== Alignment Statistics Summary ==="
column -t -s$'\t' "$OUTPUT"
echo ""
echo "Output written to: $OUTPUT"
echo "Done."
