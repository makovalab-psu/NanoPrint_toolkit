#!/bin/bash

# CONFIG.sh - Generate Snakefile from CONFIG file
#
# Usage: CONFIG.sh [-i CONFIG] [-o Snakefile]
#
# Parses the CONFIG file and generates a Snakefile with hardcoded
# wildcard values and output file lists.

set -euo pipefail

# Default values
CONFIG_FILE="CONFIG"
OUTPUT_FILE="Snakefile"

# Parse command line arguments
while getopts "i:o:h" opt; do
    case $opt in
        i) CONFIG_FILE="$OPTARG" ;;
        o) OUTPUT_FILE="$OPTARG" ;;
        h)
            echo "Usage: CONFIG.sh [-i CONFIG] [-o Snakefile]"
            echo ""
            echo "Options:"
            echo "  -i    Input CONFIG file (default: CONFIG)"
            echo "  -o    Output Snakefile (default: Snakefile)"
            echo "  -h    Show this help message"
            exit 0
            ;;
        *)
            echo "Usage: CONFIG.sh [-i CONFIG] [-o Snakefile]" >&2
            exit 1
            ;;
    esac
done

# Check CONFIG file exists
if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "Error: CONFIG file not found: $CONFIG_FILE" >&2
    exit 1
fi

# Initialize arrays
declare -a GENOMES=()
declare -a FEATURES=()
declare -a FEATURE_N_WINDOWS=()
declare -a FEATURE_WINDOW_SIZES=()
declare -a WINDOW_SIZES=()
declare -a SIG_LEVELS=()
declare -a SAMPLES=()
declare -a TREATMENTS=()
declare -a CONTROLS=()
declare -a TEMP_DIRS=()

# Function to strip file extension
strip_ext() {
    local filename="$1"
    # Remove common extensions
    filename="${filename%.gz}"
    filename="${filename%.fastq}"
    filename="${filename%.fq}"
    filename="${filename%.bam}"
    filename="${filename%.fa}"
    filename="${filename%.fasta}"
    filename="${filename%.bed}"
    echo "$filename"
}

# Parse CONFIG file
while IFS= read -r line || [[ -n "$line" ]]; do
    # Skip empty lines
    [[ -z "$line" ]] && continue

    # Skip comment lines (starting with # but not ^)
    [[ "$line" =~ ^#[^\ ] ]] && continue
    [[ "$line" =~ ^#[[:space:]] ]] && continue

    # Handle prefix lines (^X followed by tab-separated values)
    if [[ "$line" =~ ^\^ ]]; then
        prefix="${line:0:2}"
        # Get values after prefix (tab-separated)
        values="${line:2}"
        # Trim leading whitespace/tabs
        values="${values#"${values%%[![:space:]]*}"}"

        case "$prefix" in
            "^g")
                # Genome - strip extension
                value=$(echo "$values" | cut -f1)
                value=$(strip_ext "$value")
                GENOMES+=("$value")
                ;;
            "^f")
                # Feature - strip extension, get n_windows and window_size
                value=$(echo "$values" | cut -f1)
                value=$(strip_ext "$value")
                n_windows=$(echo "$values" | cut -f2 | tr -d ' ')
                window_size=$(echo "$values" | cut -f3 | tr -d ' ')
                # Use defaults if not specified
                [[ -z "$n_windows" ]] && n_windows="1000"
                [[ -z "$window_size" ]] && window_size="10"
                FEATURES+=("$value")
                FEATURE_N_WINDOWS+=("$n_windows")
                FEATURE_WINDOW_SIZES+=("$window_size")
                ;;
            "^w")
                # Window size
                value=$(echo "$values" | cut -f1 | tr -d ' ')
                WINDOW_SIZES+=("$value")
                ;;
            "^s")
                # Significance threshold
                value=$(echo "$values" | cut -f1 | tr -d ' ')
                SIG_LEVELS+=("$value")
                ;;
            "^r")
                # Relationship: sample, treatment, control
                sample=$(echo "$values" | cut -f1)
                treatment=$(echo "$values" | cut -f2)
                control=$(echo "$values" | cut -f3)
                # Strip extensions
                sample=$(strip_ext "$sample")
                treatment=$(strip_ext "$treatment")
                control=$(strip_ext "$control")
                SAMPLES+=("$sample")
                TREATMENTS+=("$treatment")
                CONTROLS+=("$control")
                ;;
            "^t")
                # Temporary directory pattern
                value=$(echo "$values" | cut -f1)
                # Trim trailing slashes
                value="${value%/}"
                TEMP_DIRS+=("$value")
                ;;
        esac
    fi
done < "$CONFIG_FILE"

# Collect all raw samples (unique treatment + control samples)
declare -a RAW_SAMPLES=()
for t in "${TREATMENTS[@]}"; do
    RAW_SAMPLES+=("$t")
done
for c in "${CONTROLS[@]}"; do
    RAW_SAMPLES+=("$c")
done
# Remove duplicates
RAW_SAMPLES=($(printf '%s\n' "${RAW_SAMPLES[@]}" | sort -u))

# ============================================================================
# Extract chromosomes from genome .fai files
# ============================================================================
# Use parallel arrays for bash 3.2 compatibility (no associative arrays)
declare -a GENOME_CHR_NAMES=()    # genome names (parallel to GENOME_CHR_VALUES)
declare -a GENOME_CHR_VALUES=()   # space-separated chromosome lists

for genome in "${GENOMES[@]}"; do
    fai_file="resources/genomes/${genome}.fa.fai"
    if [[ ! -f "$fai_file" ]]; then
        echo "Creating index for ${genome}.fa..."
        samtools faidx "resources/genomes/${genome}.fa"
    fi
    # Extract chromosome names from column 1
    chrs=$(cut -f1 "$fai_file" | tr '\n' ' ')
    GENOME_CHR_NAMES+=("$genome")
    GENOME_CHR_VALUES+=("$chrs")
done

# Helper to get chromosomes for a genome (returns via RESULT variable)
get_genome_chrs() {
    local target="$1"
    RESULT=""
    for i in "${!GENOME_CHR_NAMES[@]}"; do
        if [[ "${GENOME_CHR_NAMES[$i]}" == "$target" ]]; then
            RESULT="${GENOME_CHR_VALUES[$i]}"
            return
        fi
    done
}

# ============================================================================
# Extract chromosomes from feature BED files (filtered to genome chromosomes)
# ============================================================================
declare -a FEATURE_CHR_NAMES=()   # feature names (parallel to FEATURE_CHR_VALUES)
declare -a FEATURE_CHR_VALUES=()  # space-separated chromosome lists

for feature in "${FEATURES[@]}"; do
    bed_file="resources/features/${feature}.bed"
    if [[ ! -f "$bed_file" ]]; then
        echo "Error: Feature BED file not found: $bed_file" >&2
        exit 1
    fi

    # Get genome chromosomes (using first genome as reference)
    get_genome_chrs "${GENOMES[0]}"
    genome_chrs="$RESULT"

    # Get unique chromosomes from BED file, filter to those in genome
    bed_chrs=$(cut -f1 "$bed_file" | sort -u)
    filtered_chrs=""
    for chr in $bed_chrs; do
        # Check if chr is in genome_chrs (space-delimited search)
        if [[ " $genome_chrs " == *" $chr "* ]]; then
            filtered_chrs="$filtered_chrs $chr"
        fi
    done
    # Trim leading space
    filtered_chrs="${filtered_chrs# }"
    FEATURE_CHR_NAMES+=("$feature")
    FEATURE_CHR_VALUES+=("$filtered_chrs")
done

# Helper to get chromosomes for a feature (returns via RESULT variable)
get_feature_chrs() {
    local target="$1"
    RESULT=""
    for i in "${!FEATURE_CHR_NAMES[@]}"; do
        if [[ "${FEATURE_CHR_NAMES[$i]}" == "$target" ]]; then
            RESULT="${FEATURE_CHR_VALUES[$i]}"
            return
        fi
    done
}

# Helper function to format array as Python list
python_list() {
    local arr=("$@")
    local result="["
    local first=true
    for item in "${arr[@]}"; do
        if $first; then
            first=false
        else
            result+=", "
        fi
        result+="\"$item\""
    done
    result+="]"
    echo "$result"
}

# Generate Snakefile
cat > "$OUTPUT_FILE" << 'EOF'
# ============================================================================
# Auto-generated Snakefile by CONFIG.sh
# Do not edit manually - regenerate from CONFIG file
# ============================================================================

import os

# ============================================================================
# Configuration
# ============================================================================

EOF

# Write wildcard lists
{
    echo "# Genomes"
    echo "GENOMES = $(python_list "${GENOMES[@]}")"
    echo ""
    echo "# Features for annotation"
    echo "FEATURES = $(python_list "${FEATURES[@]}")"
    echo ""
    echo "# Feature annotation parameters: feature -> (n_windows, window_size)"
    echo "FEATURE_PARAMS = {"
    for i in "${!FEATURES[@]}"; do
        echo "    \"${FEATURES[$i]}\": (${FEATURE_N_WINDOWS[$i]}, ${FEATURE_WINDOW_SIZES[$i]}),"
    done
    echo "}"
    echo ""
    echo "# Window sizes for density calculation"
    echo "WINDOW_SIZES = $(python_list "${WINDOW_SIZES[@]}")"
    echo ""
    echo "# Significance threshold levels"
    echo "SIG_LEVELS = $(python_list "${SIG_LEVELS[@]}")"
    echo ""
    echo "# Strand directions"
    echo "STRANDS = [\"for\", \"rev\"]"
    echo ""
    echo "# All raw samples (treatment and control)"
    echo "RAW_SAMPLES = $(python_list "${RAW_SAMPLES[@]}")"
    echo ""
    echo "# Sample names (from relationships)"
    echo "SAMPLES = $(python_list "${SAMPLES[@]}")"
    echo ""
    echo "# Relationship mapping: sample -> (treatment_sample, control_sample)"
    echo "RELATIONSHIPS = {"
    for i in "${!SAMPLES[@]}"; do
        echo "    \"${SAMPLES[$i]}\": (\"${TREATMENTS[$i]}\", \"${CONTROLS[$i]}\"),"
    done
    echo "}"
    echo ""
    echo "# Chromosomes per genome (from .fai files)"
    echo "CHROMOSOMES = {"
    for i in "${!GENOME_CHR_NAMES[@]}"; do
        genome="${GENOME_CHR_NAMES[$i]}"
        chrs_array=(${GENOME_CHR_VALUES[$i]})
        echo "    \"$genome\": $(python_list "${chrs_array[@]}"),"
    done
    echo "}"
    echo ""
    echo "# Chromosomes per feature (filtered to genome chromosomes)"
    echo "FEATURE_CHROMOSOMES = {"
    for i in "${!FEATURE_CHR_NAMES[@]}"; do
        feature="${FEATURE_CHR_NAMES[$i]}"
        chrs_array=(${FEATURE_CHR_VALUES[$i]})
        echo "    \"$feature\": $(python_list "${chrs_array[@]}"),"
    done
    echo "}"
    echo ""
    echo "# Directories to mark as temporary (auto-deleted after downstream rules complete)"
    if [[ ${#TEMP_DIRS[@]} -gt 0 ]]; then
        echo "TEMP_DIRS = $(python_list "${TEMP_DIRS[@]}")"
    else
        echo "TEMP_DIRS = []"
    fi
    echo ""
    echo "# Mapping of output keys to whether they should be temporary"
    echo "TEMP_OUTPUTS = {"
    echo "    \"aligned_reads_bam\": \"data/aligned_reads\" in TEMP_DIRS,"
    echo "    \"perbase_error_by_chr\": \"data/perbase_error_by_chr\" in TEMP_DIRS,"
    echo "    \"reactivity\": \"data/reactivity\" in TEMP_DIRS,"
    echo "    \"bg\": \"data/bg\" in TEMP_DIRS,"
    echo "    \"windows\": \"data/windows\" in TEMP_DIRS,"
    echo "    \"bw\": \"data/bw\" in TEMP_DIRS,"
    echo "    \"annotations\": \"data/annotations\" in TEMP_DIRS,"
    echo "    \"annotations_merged\": \"data/annotations_merged\" in TEMP_DIRS,"
    echo "}"
    echo ""
} >> "$OUTPUT_FILE"

# Write helper functions and rules
cat >> "$OUTPUT_FILE" << 'EOF'
# ============================================================================
# Helper functions for relationship lookups
# ============================================================================

def get_treatment(sample):
    """Get treatment sample name for a relationship sample."""
    return RELATIONSHIPS[sample][0]

def get_control(sample):
    """Get control sample name for a relationship sample."""
    return RELATIONSHIPS[sample][1]

def wrap_output(key, path):
    """Wrap output path in temp() if configured as temporary."""
    if TEMP_OUTPUTS.get(key, False):
        return temp(path)
    return path

# ============================================================================
# Include rule files
# ============================================================================

include: "workflow/rules/phase1_mapping_qc.smk"
include: "workflow/rules/phase2_perbase_error.smk"
include: "workflow/rules/phase3_reactivity.smk"
include: "workflow/rules/phase4_analysis.smk"
include: "workflow/rules/phase5_annotate_features.smk"

# ============================================================================
# Target rule
# ============================================================================

rule all:
    input:
        # Phase 1: Read statistics
        expand("tables/read_stats/{raw_sample}.txt",
               raw_sample=RAW_SAMPLES),

        # Phase 1: Alignment statistics
        expand("tables/alignment_stats/{genome}/{raw_sample}.txt",
               genome=GENOMES, raw_sample=RAW_SAMPLES),

        # Phase 2: Per-base error files
        expand("data/perbase_error/{genome}/{raw_sample}_{strand}.txt.gz",
               genome=GENOMES, raw_sample=RAW_SAMPLES, strand=STRANDS),

        # Phase 4: Merged bigWig files
        expand("data/bw_merged/{genome}/significance_threshold_{sig}/{sample}_{strand}.bw",
               genome=GENOMES, sig=SIG_LEVELS, sample=SAMPLES, strand=STRANDS),

        # Phase 4: Merged density files
        expand("data/windows_merged/{genome}/window_size_{size}/significance_threshold_{sig}/{sample}_{strand}.bg",
               genome=GENOMES, size=WINDOW_SIZES, sig=SIG_LEVELS, sample=SAMPLES, strand=STRANDS),

        # Phase 5: Averaged feature annotations
        expand("data/annotations_averaged/{genome}/{feature}/{sample}_{strand}.txt.gz",
               genome=GENOMES, feature=FEATURES, sample=SAMPLES, strand=STRANDS),

EOF

echo "Generated $OUTPUT_FILE from $CONFIG_FILE"
echo ""
echo "Summary:"
echo "  Genomes:        ${#GENOMES[@]} (${GENOMES[*]})"
echo "  Features:       ${#FEATURES[@]} (${FEATURES[*]})"
echo "  Window sizes:   ${#WINDOW_SIZES[@]} (${WINDOW_SIZES[*]})"
echo "  Sig levels:     ${#SIG_LEVELS[@]} (${SIG_LEVELS[*]})"
echo "  Relationships:  ${#SAMPLES[@]}"
for i in "${!SAMPLES[@]}"; do
    echo "    ${SAMPLES[$i]}: ${TREATMENTS[$i]} vs ${CONTROLS[$i]}"
done
echo ""
echo "Chromosomes per genome:"
for i in "${!GENOME_CHR_NAMES[@]}"; do
    genome="${GENOME_CHR_NAMES[$i]}"
    chrs_array=(${GENOME_CHR_VALUES[$i]})
    echo "  $genome: ${#chrs_array[@]} chromosomes"
done
echo ""
echo "Chromosomes per feature (filtered to genome):"
for i in "${!FEATURE_CHR_NAMES[@]}"; do
    feature="${FEATURE_CHR_NAMES[$i]}"
    chrs_array=(${FEATURE_CHR_VALUES[$i]})
    echo "  $feature: ${#chrs_array[@]} chromosomes"
done
echo ""
echo "Temporary directories:"
if [[ ${#TEMP_DIRS[@]} -gt 0 ]]; then
    for dir in "${TEMP_DIRS[@]}"; do
        echo "  $dir"
    done
else
    echo "  (none)"
fi
