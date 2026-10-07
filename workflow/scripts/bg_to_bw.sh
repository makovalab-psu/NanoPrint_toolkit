#!/bin/bash

# Convert BedGraph to BigWig format
# Filters by significance threshold and converts using bedGraphToBigWig

set -euo pipefail

# Default values
SIG_LEVEL=4

# Usage function
usage() {
    cat << EOF
Usage: $(basename "$0") -i <input.bg> -o <output.bw> -g <chrom.sizes.fai> [-p sig_level] [-T tmpdir]

Convert bedGraph to bigWig format with optional significance filtering.

Required arguments:
    -i    Input bedGraph file (from react_to_bg.sh, with significance header)
    -o    Output bigWig file
    -g    Chromosome sizes file (FAI format from samtools faidx)

Optional arguments:
    -p    Significance level filter (default: 4)
              0 = return all data (no filtering)
              1 = p <= 0.05
              2 = p <= 0.01
              3 = p <= 0.001
              4 = p <= 0.0001
    -T    Temporary directory (default: same directory as output)
    -h    Show this help message

Input format: bedGraph with significance thresholds in header
    Header example:
    # p < 0.05 (black):    reactivity >= 0.002886645
    # p < 0.01 (#FF8C00):  reactivity >= 0.004327449
    ...

Output format: UCSC bigWig binary format

Example:
    $(basename "$0") -i react_chr1.bg -o react_chr1.bw -g genome.fa.fai -p 4
    $(basename "$0") -i react_chr1.bg -o react_chr1.bw -g genome.fa.fai -p 0  # all data
EOF
    exit 1
}

# Parse arguments
INPUT=""
OUTPUT=""
GENOME=""
TMP_DIR=""

while getopts "i:o:g:p:T:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        g) GENOME="$OPTARG" ;;
        p) SIG_LEVEL="$OPTARG" ;;
        T) TMP_DIR="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Validate required arguments
if [[ -z "$INPUT" || -z "$OUTPUT" || -z "$GENOME" ]]; then
    echo "Error: Missing required arguments" >&2
    usage
fi

# Check input files exist
if [[ ! -f "$INPUT" ]]; then
    echo "Error: Input file not found: $INPUT" >&2
    exit 1
fi

if [[ ! -f "$GENOME" ]]; then
    echo "Error: Genome sizes file not found: $GENOME" >&2
    exit 1
fi

# Check bedGraphToBigWig is available
if ! command -v bedGraphToBigWig &> /dev/null; then
    echo "Error: bedGraphToBigWig not found. Please install UCSC tools." >&2
    echo "Download from: http://hgdownload.soe.ucsc.edu/admin/exe/" >&2
    exit 1
fi

# Validate significance level
if [[ ! "$SIG_LEVEL" =~ ^[0-4]$ ]]; then
    echo "Error: Significance level (-p) must be 0, 1, 2, 3, or 4" >&2
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
    TMP_DIR="${OUT_DIR}/tmp_bw_$$"
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
FILTERED_BG="${TMP_DIR}/filtered.bg"
SORTED_BG="${TMP_DIR}/sorted.bg"
FILTERED_GENOME="${TMP_DIR}/filtered_genome.txt"

# Determine threshold based on significance level
if [[ "$SIG_LEVEL" -eq 0 ]]; then
    echo "Returning all data (no significance filtering)..."
    THRESHOLD="-999999"  # Will match all values
else
    # Map significance level to header pattern
    case $SIG_LEVEL in
        1) SIG_PATTERN="p < 0.05" ;;
        2) SIG_PATTERN="p < 0.01" ;;
        3) SIG_PATTERN="p < 0.001" ;;
        4) SIG_PATTERN="p < 1e-04" ;;
    esac

    echo "Parsing significance threshold from header..."
    # Extract threshold value from header
    # '|| true': when grep finds no threshold line the pipeline fails, and under
    # 'set -euo pipefail' a failing command substitution in an assignment ends the
    # script right here, before any of the messages below - the log just stops at
    # "Parsing significance threshold from header...".
    # The pattern accepts any number format R prints (0.0145786, 5e-04, 0).
    THRESHOLD=$(grep "$SIG_PATTERN" "$INPUT" | grep -oE '[0-9][0-9.e+-]*$' | head -1 || true)

    if [[ -z "$THRESHOLD" ]]; then
        DATA_COUNT=$(grep -v "^#" "$INPUT" | grep -c -v "^$" || true)
        if grep -q "No significance thresholds" "$INPUT" || [[ "$DATA_COUNT" -eq 0 ]]; then
            # react_to_bg.sh writes this header when there are no negative
            # reactivity values to build a null distribution from (or no data at
            # all). Without a null there is no threshold, so nothing can be called
            # significant: same 0-byte placeholder as "no positions passed" below.
            echo "Warning: $INPUT has no significance thresholds ($DATA_COUNT positions," >&2
            echo "  no negative reactivity values for a null distribution)." >&2
            echo "Writing an empty placeholder: $OUTPUT" >&2
            : > "$OUTPUT"
            exit 0
        fi
        echo "Error: Could not find threshold for '$SIG_PATTERN' in header" >&2
        echo "Header contents:" >&2
        head -10 "$INPUT" >&2
        exit 1
    fi

    echo "Significance level: $SIG_PATTERN"
    echo "Reactivity threshold: >= $THRESHOLD"
fi

echo "Filtering bedGraph..."
# Filter out header and values below threshold
if [[ "$SIG_LEVEL" -eq 0 ]]; then
    # No filtering - just remove header
    grep -v "^#" "$INPUT" | grep -v "^$" > "$FILTERED_BG"
else
    # Filter by threshold
    awk -v thresh="$THRESHOLD" '
        !/^#/ && !/^$/ && $4 >= thresh { print }
    ' "$INPUT" > "$FILTERED_BG"
fi

FILTERED_COUNT=$(wc -l < "$FILTERED_BG" | tr -d ' ')

if [[ "$FILTERED_COUNT" -eq 0 ]]; then
    # Not an error. On a short contig it is ordinary for no position to clear the
    # significance threshold, and a sample has next to no reads on a reference it
    # does not contain (js4031: a G4 library on the perfect reverse complement).
    # bedGraphToBigWig cannot write an empty bigWig, so leave a 0-byte placeholder,
    # as igv_coverage_bigwig does; Merge_bigwig.sh skips it.
    INPUT_COUNT=$(grep -v "^#" "$INPUT" | grep -c -v "^$" || true)
    echo "Warning: none of the $INPUT_COUNT positions in $INPUT passed the filter." >&2
    echo "Writing an empty placeholder: $OUTPUT" >&2
    : > "$OUTPUT"
    exit 0
fi

# Extract chromosomes from filtered bedGraph
echo "Extracting chromosomes from input..."
INPUT_CHRS=$(cut -f1 "$FILTERED_BG" | sort -u)

# Filter genome file to only include chromosomes from input
echo "Filtering genome file to matching chromosomes..."
awk -v chrs="$INPUT_CHRS" '
BEGIN {
    n = split(chrs, arr, "\n")
    for (i in arr) valid[arr[i]] = 1
}
$1 in valid { print $1, $2 }
' OFS='\t' "$GENOME" > "$FILTERED_GENOME"

if [[ ! -s "$FILTERED_GENOME" ]]; then
    echo "Error: No matching chromosomes found between input and genome file" >&2
    exit 1
fi

# Sort bedGraph by chromosome and position (required by bedGraphToBigWig)
echo "Sorting bedGraph..."
sort -k1,1 -k2,2n "$FILTERED_BG" > "$SORTED_BG"

# Convert to bigWig
echo "Converting to bigWig..."
bedGraphToBigWig "$SORTED_BG" "$FILTERED_GENOME" "$OUTPUT"

echo ""
echo "Done: $OUTPUT"
