#!/bin/bash

# Description: Filter alignments - keep mapped primary alignments at or above a
# mapping quality. Unmapped, secondary and supplementary records are always dropped.
# Usage: ./Filter_alignments.sh -i <input_bam> -o <output_bam> [-q <min_mapq>]
#
# -q defaults to 20. -q 0 keeps every mapped primary alignment ("unfiltered").
# Unmapped reads are excluded with a flag, not by MAPQ: they carry MAPQ 0, so at
# -q 0 a MAPQ test alone would let them through, and uncalled4 cannot align them.
# Secondary and supplementary records stay out at every setting - a supplementary
# record is hard-clipped while its move table still describes the whole read (see
# the concatemer note in claude.md).

# Function to display usage
usage() {
    echo "Usage: $0 -i <input_bam> -o <output_bam> [-q <min_mapq, default 20>]"
    exit 1
}

# Parse command line arguments
MIN_MAPQ=20

while getopts ":i:o:q:g:T:" opt; do
    case ${opt} in
        i )
            INPUT_FILE=$OPTARG
            ;;
        o )
            OUTPUT_FILE=$OPTARG
            ;;
        q )
            MIN_MAPQ=$OPTARG
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

if ! [[ "${MIN_MAPQ}" =~ ^[0-9]+$ ]]; then
    echo "Error: -q must be a non-negative integer, got '${MIN_MAPQ}'" >&2
    exit 1
fi

echo "Filtering ${INPUT_FILE} (mapped primary alignments, MAPQ >= ${MIN_MAPQ})"

# 0x904 = unmapped (0x4) + secondary (0x100) + supplementary (0x800)
samtools view -b -q "${MIN_MAPQ}" -F 0x904 "${INPUT_FILE}" > "${OUTPUT_FILE}"

echo "Filtering completed. Filtered bam file available at ${OUTPUT_FILE}"

# Cleanup function
cleanup() {
    if [[ -d "${TMP_DIR}" ]]; then
        rm -rf "${TMP_DIR}"
    fi
}
trap cleanup EXIT

