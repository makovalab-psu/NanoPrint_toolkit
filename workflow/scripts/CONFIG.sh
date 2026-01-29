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
