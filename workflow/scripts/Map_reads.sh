#!/bin/bash

# Description: Map raw reads to reference genome using minimap2 and sort with samtools
# Usage: ./Map_reads.sh -i <input_file> -o <output_file> -g <genome.fasta> -T <temp_dir>

# Function to display usage
usage() {
    echo "Usage: $0 -i <input_file> -o <output_file> -g <genome.fasta> -T <temp_dir>"
    exit 1
}

# Parse command line arguments
while getopts ":i:o:g:T:" opt; do
    case ${opt} in
        i )
            INPUT_FILE=$OPTARG
            ;;
        o )
            OUTPUT_FILE=$OPTARG
            ;;
        g )
            GENOME_FILE=$OPTARG
            ;;
        T )
            TMP_DIR=$OPTARG
            ;;
        \? )
            usage
            ;;
    esac
done

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

# Map reads and output to SAM format

if [ "${IS_BAM}" = true ]; then
    echo "Mapping BAM file ${INPUT_FILE}..."
    FASTQ_FILE="${INPUT_FILE}.bam"
    samtools fastq "${INPUT_FILE}" > "${FASTQ_FILE}"
else
    echo "Mapping FASTQ file ${INPUT_FILE}..."
    FASTQ_FILE="${INPUT_FILE}"
fi

minimap2 -a -x lr:hq "${GENOME_FILE}" "${FASTQ_FILE}" > "${TMP_DIR}/output.sam" 

# Sort the SAM file and convert to BAM
echo "Sorting the output..."
samtools sort -o "${OUTPUT_FILE}" "${TMP_DIR}/output.sam"

echo "Mapping completed. Sorted BAM file available at ${OUTPUT_FILE}."

# Cleanup function
cleanup() {
    if [[ -d "$TMP_DIR" ]]; then
        rm -rf "$TMP_DIR"
    fi
}
trap cleanup EXIT

