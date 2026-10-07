#!/bin/bash
set -euo pipefail

# Description: Map raw reads to reference genome using minimap2 and sort with samtools
# Usage: ./Map_reads.sh -i <input_file> -o <output_file> -g <genome.fasta> -T <temp_dir> -t <threads>
# When input is a BAM (e.g. dorado basecalled), MM/ML/MN/mv/ts/pi/sp/ns/fn/BC/RG/qs tags are
# preserved via samtools fastq -T and minimap2 -y so downstream tools can use them.
# mv/ts/pi/sp/ns are the move table Uncalled4 needs; fn (source pod5 basename) is needed by
# nanoprint preprocess's per-pod5 batching split; MM/ML are the base modification calls and
# MN the sequence length they were called against (lets modkit validate them); BC is the
# dorado barcode classification that nanoprint demux splits on; RG is the read group
# (run + model + barcode); qs is the mean read Q-score.
#
# THIS LIST IS AN ALLOWLIST: any tag not named here is silently dropped at alignment. That
# has already caused two real data-loss bugs — `fn` (js4017, every pod5 split came back
# empty) and MM/ML (js4014, a whole WGS BAM reached modkit with no methylation calls in it
# and every record failed). If a new consumer needs another dorado tag, add it HERE and
# grep for the tag-list comments that mirror this one.
#
# minimap2 also writes its own header, which has no @RG lines, so per-read RG tags would
# point at read groups the header does not define. For BAM input the input's @RG header
# lines are fed into the SAM stream ahead of minimap2's output (header lines may come in any
# order as long as they precede the records), so the run/flowcell/model metadata survives.

# Function to display usage
usage() {
    echo "Usage: $0 -i <input_file> -o <output_file> -g <genome.fasta> -T <temp_dir> -t <threads>"
    exit 1
}

# Parse command line arguments
INPUT_FILE=""
OUTPUT_FILE=""
GENOME_FILE=""
TMP_DIR=""
THREADS=""

while getopts ":i:o:g:T:t:" opt; do
    case ${opt} in
        i ) INPUT_FILE=$OPTARG ;;
        o ) OUTPUT_FILE=$OPTARG ;;
        g ) GENOME_FILE=$OPTARG ;;
        T ) TMP_DIR=$OPTARG ;;
        t ) THREADS=$OPTARG ;;
        \? ) usage ;;
    esac
done

THREADS="${THREADS:-1}"

# Check if all mandatory arguments are provided
if [ -z "${INPUT_FILE}" ] || [ -z "${OUTPUT_FILE}" ] || [ -z "${GENOME_FILE}" ]; then
    usage
fi

# Determine input file type
if [[ "${INPUT_FILE}" == *.bam ]]; then
    IS_BAM=true
elif [[ "${INPUT_FILE}" == *.fastq.gz ]]; then
    IS_BAM=false
else
    echo "Error: Input file must be a .bam or .fastq.gz file."
    exit 1
fi

# Get output directory
OUT_DIR=$(dirname "$OUTPUT_FILE")
if [[ -z "$OUT_DIR" || "$OUT_DIR" == "." ]]; then
    OUT_DIR="$(pwd)"
fi
mkdir -p "$OUT_DIR"

# Set temp directory (default: in output directory)
if [[ -z "$TMP_DIR" ]]; then
    TMP_DIR="${OUT_DIR}/tmp_mapping_$$"
fi
mkdir -p "$TMP_DIR"

cleanup() { [[ -d "$TMP_DIR" ]] && rm -rf "$TMP_DIR"; }
trap cleanup EXIT

# Map reads and output to BAM format
# BAM input: pipe samtools fastq (preserving dorado tags) → minimap2 -y (propagate tags) → sort
# FASTQ input: pipe minimap2 → sort directly
echo "Sorting the output..."
if [ "${IS_BAM}" = true ]; then
    echo "Mapping BAM file ${INPUT_FILE} (preserving MM/ML mod calls, mv/ts/fn, BC/RG tags)..."
    RG_HEADER="${TMP_DIR}/rg_header.sam"
    samtools view -H "${INPUT_FILE}" | { grep '^@RG' || true; } > "${RG_HEADER}"
    echo "Carrying $(wc -l < "${RG_HEADER}" | tr -d ' ') @RG header line(s) into the aligned BAM"
    {
        cat "${RG_HEADER}"
        samtools fastq -T "MM,ML,MN,mv,ts,pi,sp,ns,fn,BC,RG,qs" "${INPUT_FILE}" \
            | minimap2 -y -a -x lr:hq -t "${THREADS}" "${GENOME_FILE}" -
    } | samtools sort -@ "${THREADS}" -T "${TMP_DIR}/sort_tmp" -o "${OUTPUT_FILE}"
else
    echo "Mapping FASTQ file ${INPUT_FILE}..."
    minimap2 -a -x lr:hq -t "${THREADS}" "${GENOME_FILE}" "${INPUT_FILE}" \
        | samtools sort -@ "${THREADS}" -T "${TMP_DIR}/sort_tmp" -o "${OUTPUT_FILE}"
fi

echo "Mapping completed. Sorted BAM file available at ${OUTPUT_FILE}."

