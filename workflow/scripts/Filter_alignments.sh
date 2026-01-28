#!/bin/bash

# Description: Map raw reads to reference genome using minimap2 and sort with samtools
# Usage: ./Filter_alignments.sh -i <input_bam> -o <output_file> -g <genome.fasta> -T <temp_dir>

# Function to display usage
usage() {
    echo "Usage: $0 -i <input_bam> -o <output_bam>"
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
        \? )
            usage
            ;;
    esac
done

# Check if all mandatory arguments are provided
if [ -z "${INPUT_FILE}" ] || [ -z "${OUTPUT_FILE}" ] ; then
    usage
fi

echo "Filtering ${INPUT_FILE}"

samtools view -b -q 20 -F 0x100 -F 0x800  "${INPUT_FILE}" > "${OUTPUT_FILE}"

echo "Filtering completed. Filtered bam file available at ${OUTPUT_FILE}"

# Cleanup function
cleanup() {
    if [[ -d "${TMP_DIR}" ]]; then
        rm -rf "${TMP_DIR}"
    fi
}
trap cleanup EXIT

