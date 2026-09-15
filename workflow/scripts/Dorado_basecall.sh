#!/bin/bash

# Dorado_basecall.sh - Basecall pod5 files using Dorado
# Produces a BAM with --emit-moves so Uncalled4 can use the move table.
#
# The model string is passed to dorado verbatim, so it may carry dorado's inline
# modification syntax (e.g. sup,5mCG_5hmCG) or a full combined model name. Callers
# default to a model with CpG
# 5mC/5hmC calling enabled — see DEFAULT_DORADO_MODEL in bin/nanoprint and
# DORADO_MODEL in CONFIG.sh. A model with no modification suffix produces a BAM
# with no MM/ML tags, which cannot be recovered without basecalling again.
#
# -k <kit> classifies barcodes during basecalling (dorado --kit-name) and adds --no-trim.
# Trimming is off so the read sequence stays consistent with the move table and raw signal
# that Uncalled4 aligns; minimap2 soft-clips the barcode and adapter instead. All barcodes
# go into the one output BAM, each read tagged BC:Z:<kit>_barcodeNN (unclassified reads
# carry no BC tag). Tested with SQK-RBK114-24 on dorado 1.3.2 (js4022): mv/ts survive.
#
# Usage: Dorado_basecall.sh -i <pod5_dir_or_file> -o <output.bam> -m <model> [-t <threads>] [-k <kit>]

set -euo pipefail

usage() {
    cat << EOF
Usage: $(basename "$0") -i <pod5_dir_or_file> -o <output.bam> -m <model> [-t <threads>] [-k <kit>]

Basecall Oxford Nanopore pod5 files using Dorado.
Emits move tables (--emit-moves) required by Uncalled4 signal alignment.

Required arguments:
    -i    Input: pod5 directory or single pod5 file
    -o    Output BAM file (unsorted; sorted downstream by map_reads)
    -m    Dorado model. Accepts dorado's inline modification syntax, e.g.
          dna_r10.4.1_e8.2_400bps_sup@v5.2.0_5mCG_5hmCG@v2  (exact model, CpG calls; note the underscore)
          sup,5mCG_5hmCG                                    (auto-select + CpG calls)
          sup / hac / fast                                (no modification calls)

Optional arguments:
    -t    Number of threads (default: 1; GPU usage controlled by dorado itself)
    -k    Barcoding kit (e.g. SQK-RBK114-24). Adds --kit-name <kit> --no-trim, so reads
          are classified into BC:Z: tags in the single output BAM, untrimmed.
    -h    Show this help message

Notes:
    - Requires 'dorado' in PATH (https://github.com/nanoporetech/dorado)
    - GPU is used automatically if available; set CUDA_VISIBLE_DEVICES to control
    - Output BAM is coordinate-unsorted (read order); map_reads sorts it after alignment

Example:
    $(basename "$0") -i raw_data/Sample01/ -o data/basecalled/Sample01.bam \\
        -m dna_r10.4.1_e8.2_400bps_sup@v5.2.0_5mCG_5hmCG@v2 -t 4
EOF
    exit 1
}

INPUT=""
OUTPUT=""
MODEL=""
THREADS=1
KIT=""

while getopts "i:o:m:t:k:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        m) MODEL="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        k) KIT="$OPTARG" ;;
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
echo "Kit:     ${KIT:-none (no barcode classification)}"
echo ""

DORADO_FLAGS=(--emit-moves)
[[ -d "$INPUT" ]] && DORADO_FLAGS+=(--recursive)
[[ -n "$KIT" ]] && DORADO_FLAGS+=(--kit-name "$KIT" --no-trim)

dorado basecaller \
    "$MODEL" \
    "$INPUT" \
    "${DORADO_FLAGS[@]}" \
    > "$OUTPUT"

echo ""
echo "Done. Output: $OUTPUT"
echo "Read count: $(samtools view -c "$OUTPUT")"
