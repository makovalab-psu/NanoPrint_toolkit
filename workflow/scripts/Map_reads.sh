#!/bin/bash
set -euo pipefail

# Description: Map raw reads to reference genome using minimap2 and sort with samtools
# Usage: ./Map_reads.sh -i <input_file> -o <output_file> -g <genome.fasta> -T <temp_dir> -t <threads>
# When input is a BAM (e.g. dorado basecalled), mv/ts/pi/sp/ns tags are preserved via
# samtools fastq -T and minimap2 -y so Uncalled4 can use them downstream.

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
    echo "Mapping BAM file ${INPUT_FILE} (preserving mv/ts tags for Uncalled4)..."
    samtools fastq -T "mv,ts,pi,sp,ns" "${INPUT_FILE}" \
        | minimap2 -y -a -x lr:hq -t "${THREADS}" "${GENOME_FILE}" - \
        | samtools sort -@ "${THREADS}" -T "${TMP_DIR}/sort_tmp" -o "${OUTPUT_FILE}"
else
    echo "Mapping FASTQ file ${INPUT_FILE}..."
    minimap2 -a -x lr:hq -t "${THREADS}" "${GENOME_FILE}" "${INPUT_FILE}" \
        | samtools sort -@ "${THREADS}" -T "${TMP_DIR}/sort_tmp" -o "${OUTPUT_FILE}"
fi

echo "Mapping completed. Sorted BAM file available at ${OUTPUT_FILE}."

