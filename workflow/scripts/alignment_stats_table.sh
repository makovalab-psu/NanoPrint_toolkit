#!/bin/bash

# Alignment Stats Table: Combine per-sample alignment statistics into a single CSV table
# Pivots long-format per-file stats into wide-format with two rows per sample:
# one row for raw (not filtered) alignments, one row for filtered alignments.

set -euo pipefail

usage() {
    cat << EOF
Usage: $(basename "$0") -o <output.csv> <input1.txt> [input2.txt ...]

Combine per-sample alignment statistics into a single CSV table.
Each input file (output of Alignment_stats.sh) produces two rows:
one for raw alignments and one for filtered alignments.

Required arguments:
    -o    Output CSV file

Positional arguments:
    One or more alignment_stats txt files (tab-delimited, output of Alignment_stats.sh)

Input format (tab-delimited):
    Sample  Statistic  Raw_alignment  Filtered_alignment

Output format (CSV):
    Sample, Filter_status, Total_sequences, Total_length, Bases_mapped,
    Bases_mapped_cigar, Mismatches, Error_rate, Average_length, Average_quality,
    Primary_alignments, Secondary_alignments, Supplementary_alignments

Example:
    $(basename "$0") -o tables/alignment_stats_table.csv \\
        tables/alignment_stats/genome1/sample1.txt \\
        tables/alignment_stats/genome1/sample2.txt
EOF
    exit 1
}

OUTPUT=""

while getopts "o:h" opt; do
    case $opt in
        o) OUTPUT="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done
shift $((OPTIND - 1))

INPUTS=("$@")

if [[ -z "$OUTPUT" || ${#INPUTS[@]} -eq 0 ]]; then
    echo "Error: Missing required arguments" >&2
    usage
fi

mkdir -p "$(dirname "$OUTPUT")"

echo "Combining ${#INPUTS[@]} alignment stats files..."

# Write CSV header
echo "Sample,Filter_status,Total_sequences,Total_length,Bases_mapped,Bases_mapped_cigar,Mismatches,Error_rate,Average_length,Average_quality,Primary_alignments,Secondary_alignments,Supplementary_alignments" > "$OUTPUT"

# Process each input file
for f in "${INPUTS[@]}"; do
    if [[ ! -f "$f" ]]; then
        echo "Warning: File not found: $f" >&2
        continue
    fi

    echo "  Processing: $f"

    # Pivot long-format stats into two wide-format rows (raw and filtered).
    # Input columns (tab-delimited): Sample, Statistic, Raw_alignment, Filtered_alignment
    # Output: two CSV rows per file — "Not_filtered" using col 3, "Filtered" using col 4.
    tail -n +2 "$f" | awk '
    BEGIN {
        FS = "\t"
        OFS = ","
        # Ordered list of stat names to extract (must match Alignment_stats.sh output)
        n = 11
        want[1]  = "Total sequences"
        want[2]  = "Total length"
        want[3]  = "Bases mapped"
        want[4]  = "Bases mapped (cigar)"
        want[5]  = "Mismatches"
        want[6]  = "Error rate"
        want[7]  = "Average length"
        want[8]  = "Average quality"
        want[9]  = "Primary alignments"
        want[10] = "Secondary alignments"
        want[11] = "Supplementary alignments"
    }
    {
        sample   = $1
        stat     = $2
        raw[stat]      = $3
        filtered[stat] = $4
    }
    END {
        # Row 1: not filtered (raw alignment values)
        printf "%s,%s", sample, "Not_filtered"
        for (i = 1; i <= n; i++) {
            printf ",%s", raw[want[i]]
        }
        printf "\n"

        # Row 2: filtered alignment values
        printf "%s,%s", sample, "Filtered"
        for (i = 1; i <= n; i++) {
            printf ",%s", filtered[want[i]]
        }
        printf "\n"
    }
    ' >> "$OUTPUT"
done

echo "Done: $OUTPUT"
