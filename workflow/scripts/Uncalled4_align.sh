#!/bin/bash

# Uncalled4_align.sh - Align raw nanopore signals to the pore model using Uncalled4
# Takes a reference-aligned BAM with dorado move tags (mv,ts,pi,sp,ns) and pod5 files.
# The BAM must be aligned AND retain the move tags from dorado --emit-moves.
# Map_reads.sh preserves these tags via samtools fastq -T "mv,ts,pi,sp,ns" + minimap2 -y.
# Produces a BAM with signal-level DTW alignment data.
#
# --min-aln-length 50: uncalled4 default is ~200 bp (for WGS). Short synthetic targets
# (G4 oligos ~86-89 bp) would all fail the default threshold. 50 bp works for both.
#
# Usage: Uncalled4_align.sh -i <filtered_alignments.bam> -p <pod5_path> -g <genome.fa> -o <out.bam> [-t <threads>]
#
# Pod5 path may be a single .pod5 file or a directory. For directories, the script
# uses find to enumerate all .pod5 files recursively into a temp list and passes that
# to uncalled4 --reads (uncalled4 does not recurse into subdirectories by default).
# Command syntax: uncalled4 align --bam-in <bam> --ref <fa> --reads <pod5> -o <out.bam> -p <n>
# Parallelism via -p (official docs: default is 1 process).

set -euo pipefail

usage() {
    cat << EOF
Usage: $(basename "$0") -i <filtered_alignments.bam> -p <pod5_path> -g <genome.fa> -o <out.bam> [-t <threads>]

Align raw nanopore signals to the pore model reference using Uncalled4 (BAM output).
Input BAM must be a reference-aligned BAM that retains dorado move tags (mv, ts, pi, sp, ns).
Map_reads.sh produces this via: samtools fastq -T "mv,ts,pi,sp,ns" | minimap2 -y -a ...
The tags survive through filter_alignments (samtools view preserves all BAM tags by default).

Required arguments:
    -i    Reference-aligned BAM with move tags (data/filtered_alignments/{genome}/{sample}.bam)
    -p    Absolute path to pod5 directory or single pod5 file
    -g    Reference genome FASTA
    -o    Output Uncalled4 BAM with DTW signal alignment

Optional arguments:
    -t    Number of parallel processes for Uncalled4 (default: 8)
    -h    Show this help message

Example:
    $(basename "$0") \\
        -i data/filtered_alignments/genome/Sample01.bam \\
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
    --bam-in         "$INPUT_BAM" \
    --ref            "$GENOME" \
    --reads          "$POD5_INPUT" \
    -p               "$THREADS" \
    --min-aln-length 50 \
    -o               "$OUTPUT" || true

if [[ ! -f "$OUTPUT" ]]; then
    echo "Error: Uncalled4 produced no BAM output" >&2
    exit 1
fi

READ_COUNT=$(samtools view -c "$OUTPUT")
if [[ "$READ_COUNT" -eq 0 ]]; then
    echo "Warning: Uncalled4 BAM is empty (0 reads)."
    echo "This is expected when the sample has no reads mapping to this genome."
fi

# uncalled4 writes reads in processing order, not coordinate order — sort before indexing.
echo ""
echo "Sorting output BAM..."
SORTED_TMP="${TMP_DIR}/sorted.bam"
samtools sort -@ "$THREADS" -T "${TMP_DIR}/sort_tmp2" -o "$SORTED_TMP" "$OUTPUT"
mv "$SORTED_TMP" "$OUTPUT"

echo "Indexing output BAM..."
samtools index "$OUTPUT"

echo "Done. Output: $OUTPUT"
echo "Read count: $READ_COUNT"
