#!/bin/bash

# Uncalled4_align.sh - Align raw nanopore signals to the pore model using Uncalled4
# Takes the dorado basecalled BAM (with move table, --emit-moves) and the original pod5 files.
# DO NOT pass a minimap2-aligned BAM — minimap2 strips the dorado mv tag that uncalled4 needs.
# uncalled4 performs its own internal alignment to the reference genome.
# Produces a BAM with signal-level DTW alignment data.
#
# Usage: Uncalled4_align.sh -i <basecalled.bam> -p <pod5_path> -g <genome.fa> -o <out.bam> [-t <threads>]
#
# Pod5 path may be a single .pod5 file or a directory. For directories, the script
# uses find to enumerate all .pod5 files recursively into a temp list and passes that
# to uncalled4 --reads (uncalled4 does not recurse into subdirectories by default).
# Command syntax: uncalled4 align --bam-in <bam> --ref <fa> --reads <pod5> -o <out.bam> -p <n>
# Parallelism via -p (official docs: default is 1 process).

set -euo pipefail

usage() {
    cat << EOF
Usage: $(basename "$0") -i <basecalled.bam> -p <pod5_path> -g <genome.fa> -o <out.bam> [-t <threads>]

Align raw nanopore signals to the pore model reference using Uncalled4 (BAM output).
Input BAM must be the dorado basecalled BAM (data/basecalled/{sample}.bam) with the
move table (--emit-moves) intact. DO NOT use a minimap2-aligned BAM — minimap2 strips
the mv tag, causing "moves missing" for all reads. uncalled4 aligns to the reference
internally.

Required arguments:
    -i    Input dorado basecalled BAM with move table (data/basecalled/{sample}.bam)
    -p    Absolute path to pod5 directory or single pod5 file
    -g    Reference genome FASTA
    -o    Output Uncalled4 BAM with DTW signal alignment

Optional arguments:
    -t    Number of parallel processes for Uncalled4 (default: 8)
    -h    Show this help message

Example:
    $(basename "$0") \\
        -i data/basecalled/Sample01.bam \\
        -p /absolute/path/to/pod5/run01/ \\
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

TMP_DIR="${OUT_DIR}/tmp_$$"
mkdir -p "$TMP_DIR"
cleanup() { [[ -d "$TMP_DIR" ]] && rm -rf "$TMP_DIR"; }
trap cleanup EXIT

echo "=== Uncalled4 Signal Alignment (BAM) ==="
echo "Input BAM:   $INPUT_BAM"
echo "Pod5 path:   $POD5_DIR"
echo "Genome:      $GENOME"
echo "Output:      $OUTPUT"
echo "Processes:   $THREADS"
echo ""

POD5_INPUT="$POD5_DIR"
if [[ -d "$POD5_DIR" ]]; then
    POD5_LIST="$TMP_DIR/pod5_files.txt"
    find "$POD5_DIR" -name "*.pod5" -type f | sort > "$POD5_LIST"
    if [[ ! -s "$POD5_LIST" ]]; then
        echo "Error: no .pod5 files found under $POD5_DIR" >&2
        exit 1
    fi
    echo "Found $(wc -l < "$POD5_LIST") pod5 file(s) (recursive search)"
    POD5_INPUT="$POD5_LIST"
fi

# uncalled4 returns non-zero when any reads fail DTW (even if most succeed).
# Use || true and verify the output is non-empty (pattern from js4007/Snakefile).
uncalled4 align \
    --bam-in "$INPUT_BAM" \
    --ref    "$GENOME" \
    --reads  "$POD5_INPUT" \
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
