#!/bin/bash

# Uncalled4_align.sh - Align raw nanopore signals to the pore model using Uncalled4
# Takes a sequence-aligned BAM (filtered_alignments) and the original pod5 files.
# Produces a BAM with signal-level DTW alignment data.
#
# Usage: Uncalled4_align.sh -i <filtered.bam> -p <pod5_dir> -g <genome.fa> -o <out.bam> [-t <threads>]
#
# Command syntax confirmed from js4004 Snakefile:
#   uncalled4 align --bam-in <bam> --ref <fa> --reads <pod5> -o <out.bam>
# Parallelism via -p (official docs: default is 1 process).

set -euo pipefail

usage() {
    cat << EOF
Usage: $(basename "$0") -i <filtered.bam> -p <pod5_dir> -g <genome.fa> -o <out.bam> [-t <threads>]

Align raw nanopore signals to the pore model reference using Uncalled4 (BAM output).
Input BAM must have sequence-level alignments (from minimap2) and a move table
(--emit-moves from dorado) so Uncalled4 can trace each read's signal back to
the raw pod5 data.

Required arguments:
    -i    Input sequence-aligned BAM (filtered_alignments/{genome}/{sample}.bam)
    -p    Pod5 file directory (raw_data/{sample}/) or single pod5 file
    -g    Reference genome FASTA
    -o    Output Uncalled4 BAM with DTW signal alignment

Optional arguments:
    -t    Number of parallel processes for Uncalled4 (default: 8)
    -h    Show this help message

Example:
    $(basename "$0") \\
        -i data/filtered_alignments/genome/Sample01.bam \\
        -p raw_data/Sample01/ \\
        -g resources/genomes/genome.fa \\
        -o data/uncalled4/genome/Sample01.bam \\
        -t 8
EOF
    exit 1
}

INPUT_BAM=""
POD5_DIR=""
GENOME=""
OUTPUT=""
THREADS=8

while getopts "i:p:g:o:t:h" opt; do
    case $opt in
        i) INPUT_BAM="$OPTARG" ;;
        p) POD5_DIR="$OPTARG" ;;
        g) GENOME="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$INPUT_BAM" || -z "$POD5_DIR" || -z "$GENOME" || -z "$OUTPUT" ]]; then
    echo "Error: Missing required arguments" >&2
    usage
fi

if [[ ! -f "$INPUT_BAM" ]]; then
    echo "Error: Input BAM not found: $INPUT_BAM" >&2
    exit 1
fi

if [[ ! -e "$POD5_DIR" ]]; then
    echo "Error: Pod5 path not found: $POD5_DIR" >&2
    exit 1
fi

if [[ ! -f "$GENOME" ]]; then
    echo "Error: Genome FASTA not found: $GENOME" >&2
    exit 1
fi

if ! command -v uncalled4 &> /dev/null; then
    echo "Error: uncalled4 not found in PATH. Install with: pip install uncalled4" >&2
    exit 1
fi

OUT_DIR=$(dirname "$OUTPUT")
mkdir -p "$OUT_DIR"

echo "=== Uncalled4 Signal Alignment (BAM) ==="
echo "Input BAM:   $INPUT_BAM"
echo "Pod5 path:   $POD5_DIR"
echo "Genome:      $GENOME"
echo "Output:      $OUTPUT"
echo "Processes:   $THREADS"
echo ""

# uncalled4 returns non-zero when any reads fail DTW (even if most succeed).
# Use || true and verify the output is non-empty (pattern from js4007/Snakefile).
uncalled4 align \
    --bam-in "$INPUT_BAM" \
    --ref    "$GENOME" \
    --reads  "$POD5_DIR" \
    -p       "$THREADS" \
    -o       "$OUTPUT" || true

if [[ ! -f "$OUTPUT" ]]; then
    echo "Error: Uncalled4 produced no BAM output" >&2
    exit 1
fi

READ_COUNT=$(samtools view -c "$OUTPUT")
if [[ "$READ_COUNT" -eq 0 ]]; then
    echo "Error: Uncalled4 BAM is empty (0 reads)" >&2
    exit 1
fi

echo ""
echo "Indexing output BAM..."
samtools index "$OUTPUT"

echo "Done. Output: $OUTPUT"
echo "Read count: $READ_COUNT"
