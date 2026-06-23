#!/bin/bash

# Dorado_basecall.sh - Basecall pod5 files using Dorado
# Produces a BAM with --emit-moves so Uncalled4 can use the move table.
#
# Usage: Dorado_basecall.sh -i <pod5_dir_or_file> -o <output.bam> -m <model> [-t <threads>]

set -euo pipefail

usage() {
    cat << EOF
Usage: $(basename "$0") -i <pod5_dir_or_file> -o <output.bam> -m <model> [-t <threads>]

Basecall Oxford Nanopore pod5 files using Dorado.
Emits move tables (--emit-moves) required by Uncalled4 signal alignment.

Required arguments:
    -i    Input: pod5 directory or single pod5 file
    -o    Output BAM file (unsorted; sorted downstream by map_reads)
    -m    Dorado model (e.g. dna_r10.4.1_e8.2_400bps_sup@v4.3.0, or 'sup', 'hac', 'fast')

Optional arguments:
    -t    Number of threads (default: 1; GPU usage controlled by dorado itself)
    -h    Show this help message

Notes:
    - Requires 'dorado' in PATH (https://github.com/nanoporetech/dorado)
    - GPU is used automatically if available; set CUDA_VISIBLE_DEVICES to control
    - Output BAM is coordinate-unsorted (read order); map_reads sorts it after alignment

Example:
    $(basename "$0") -i raw_data/Sample01/ -o data/basecalled/Sample01.bam -m sup -t 4
EOF
    exit 1
}

INPUT=""
OUTPUT=""
MODEL=""
THREADS=1

while getopts "i:o:m:t:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        m) MODEL="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$INPUT" || -z "$OUTPUT" || -z "$MODEL" ]]; then
    echo "Error: Missing required arguments" >&2
    usage
fi

if [[ ! -e "$INPUT" ]]; then
    echo "Error: Input not found: $INPUT" >&2
    exit 1
fi

if ! command -v dorado &> /dev/null; then
    echo "Error: dorado not found in PATH. Install from https://github.com/nanoporetech/dorado" >&2
    exit 1
fi

OUT_DIR=$(dirname "$OUTPUT")
mkdir -p "$OUT_DIR"

echo "=== Dorado Basecalling ==="
echo "Input:   $INPUT"
echo "Model:   $MODEL"
echo "Output:  $OUTPUT"
echo "Threads: $THREADS"
echo ""

dorado basecaller \
    "$MODEL" \
    "$INPUT" \
    --emit-moves \
    > "$OUTPUT"

echo ""
echo "Done. Output: $OUTPUT"
echo "Read count: $(samtools view -c "$OUTPUT")"
