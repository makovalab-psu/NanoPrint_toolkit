#!/bin/bash

# Description: Raise the MAPQ floor of an Uncalled4 BAM that was filtered at a lower
# MAPQ (a ^u input), keeping every tag - the DTW tags included. Used by the
# modification-calling sweep to put a ^u sample into sets at a higher MAPQ without
# re-aligning signal: Uncalled4 aligns each read on its own, so the reads that pass
# MAPQ >= N carry the same DTW whichever other reads were aligned with them.
# Usage: ./Subset_uncalled4_mapq.sh -i <uncalled4.bam> -o <output.bam> -q <min_mapq> [-t <threads>]
#
# The input must be indexed and coordinate-sorted (a ^u BAM is both); the output is
# too. Secondary, supplementary and unmapped records are dropped again, which is a
# no-op on a filtered BAM.
#
# A filter that removes nothing is reported, not failed: a genome always has some
# low-MAPQ reads (rRNA operons, IS elements), so keeping them all means either the
# data has none or the MAPQ field did not survive into this BAM - check which.

set -euo pipefail

usage() {
    echo "Usage: $0 -i <uncalled4.bam> -o <output.bam> -q <min_mapq> [-t <threads, default 1>]"
    exit 1
}

THREADS=1
while getopts ":i:o:q:t:h" opt; do
    case ${opt} in
        i ) INPUT_FILE=$OPTARG ;;
        o ) OUTPUT_FILE=$OPTARG ;;
        q ) MIN_MAPQ=$OPTARG ;;
        t ) THREADS=$OPTARG ;;
        h | \? ) usage ;;
    esac
done

if [[ -z "${INPUT_FILE:-}" || -z "${OUTPUT_FILE:-}" || -z "${MIN_MAPQ:-}" ]]; then
    usage
fi
if ! [[ "$MIN_MAPQ" =~ ^[0-9]+$ ]]; then
    echo "Error: -q must be a non-negative integer, got '$MIN_MAPQ'." >&2
    exit 1
fi
[[ -f "$INPUT_FILE" ]] || { echo "Error: input BAM not found: $INPUT_FILE" >&2; exit 1; }
[[ -f "${INPUT_FILE}.bai" ]] || { echo "Error: input BAM has no index: ${INPUT_FILE}.bai" >&2; exit 1; }

mkdir -p "$(dirname "$OUTPUT_FILE")"

# Counts from the index: instant, and cannot SIGPIPE (see the igv_coverage_bigwig note).
count_reads() {
    samtools idxstats "$1" | awk '{n += $3 + $4} END {print n + 0}'
}

N_IN=$(count_reads "$INPUT_FILE")
echo "Input:  $INPUT_FILE ($N_IN records)"
echo "Filter: MAPQ >= $MIN_MAPQ, no unmapped/secondary/supplementary"

# Written beside the output and moved into place only once indexed, so a killed job
# never leaves a BAM that looks finished.
TMP_OUT="${OUTPUT_FILE}.tmp.bam"
samtools view -@ "$THREADS" -b -h -q "$MIN_MAPQ" -F 0x904 -o "$TMP_OUT" "$INPUT_FILE"
samtools index -@ "$THREADS" "$TMP_OUT"
N_OUT=$(count_reads "$TMP_OUT")

echo "Output: $OUTPUT_FILE ($N_OUT records)"
if [[ "$N_IN" -gt 0 ]]; then
    echo "Kept:   $(awk -v a="$N_OUT" -v b="$N_IN" 'BEGIN {printf "%.2f", 100 * a / b}')%"
fi
if [[ "$N_IN" -gt 0 && "$N_OUT" -eq 0 ]]; then
    echo "Error: no record reached MAPQ $MIN_MAPQ. Either the MAPQ field was not kept in" >&2
    echo "       this Uncalled4 BAM, or the sample is entirely multimapped. Check with" >&2
    echo "       samtools view $INPUT_FILE | cut -f5 | sort -n | uniq -c" >&2
    rm -f "$TMP_OUT" "${TMP_OUT}.bai"
    exit 1
fi
mv "$TMP_OUT" "$OUTPUT_FILE"
mv "${TMP_OUT}.bai" "${OUTPUT_FILE}.bai"
if [[ "$MIN_MAPQ" -gt 0 && "$N_OUT" -eq "$N_IN" && "$N_IN" -gt 0 ]]; then
    echo "WARNING: the filter removed nothing. On a genome some reads are always below"
    echo "         MAPQ $MIN_MAPQ, so check the MAPQ column survived Uncalled4:"
    echo "         samtools view $INPUT_FILE | cut -f5 | sort -n | uniq -c"
fi
