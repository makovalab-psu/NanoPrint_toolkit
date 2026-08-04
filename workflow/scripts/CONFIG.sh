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
declare -a MEAN_WINDOW_SIZES=()
declare -a SIG_LEVELS=()
declare -a SAMPLES=()
declare -a TREATMENTS=()
declare -a CONTROLS=()
declare -a TREATMENT_PATHS=()
declare -a CONTROL_PATHS=()
declare -a TEMP_DIRS=()
IGV_BAM="false"
IGV_BIGWIG="false"
# Pinned exact model with CpG 5mC/5hmC calling enabled. The pin keeps the
# basecaller a fixed quantity across runs (the 'sup' shorthand resolves against
# each pod5's chemistry metadata, so it can drift between runs); the
# ,5mCG_5hmCG suffix is what makes dorado emit MM/ML tags at all, without which
# no methylation analysis is possible and the expensive basecalling step has to
# be repeated. A pod5 from a different chemistry will now be mis-called or
# rejected rather than auto-matched — override with '^dorado-model sup,5mCG_5hmCG'
# for auto-selection, or '^dorado-model sup' for the old mod-free behaviour.
# KEEP IN SYNC with DEFAULT_DORADO_MODEL in bin/nanoprint.
DORADO_MODEL="dna_r10.4.1_e8.2_400bps_sup@v5.2.0,5mCG_5hmCG"

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
        prefix="${line%%[[:space:]]*}"
        # Get values after prefix (tab-separated)
        values="${line#"$prefix"}"
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
                # Window size for density calculation
                value=$(echo "$values" | cut -f1 | tr -d ' ')
                WINDOW_SIZES+=("$value")
                ;;
            "^a")
                # Window size for mean reactivity
                value=$(echo "$values" | cut -f1 | tr -d ' ')
                MEAN_WINDOW_SIZES+=("$value")
                ;;
            "^s")
                # Significance threshold
                value=$(echo "$values" | cut -f1 | tr -d ' ')
                SIG_LEVELS+=("$value")
                ;;
            "^r")
                # Relationship: sample, treatment path, control path
                # treatment/control fields are absolute or relative paths (file or directory).
                # raw_sample name = basename of the path minus any recognized extension.
                sample=$(echo "$values" | cut -f1)
                treatment_raw=$(echo "$values" | cut -f2)
                control_raw=$(echo "$values" | cut -f3)
                # Strip trailing slash, then derive name from basename
                treatment=$(strip_ext "$(basename "${treatment_raw%/}")")
                control=$(strip_ext "$(basename "${control_raw%/}")")
                SAMPLES+=("$sample")
                TREATMENTS+=("$treatment")
                CONTROLS+=("$control")
                TREATMENT_PATHS+=("${treatment_raw%/}")
                CONTROL_PATHS+=("${control_raw%/}")
                ;;
            "^t")
                # Temporary directory pattern
                value=$(echo "$values" | cut -f1)
                # Trim trailing slashes
                value="${value%/}"
                TEMP_DIRS+=("$value")
                ;;
            "^igv-bam")
                IGV_BAM="true"
                ;;
            "^igv-bigwig")
                IGV_BIGWIG="true"
                ;;
            "^dorado-model")
                DORADO_MODEL="$(echo "$values" | cut -f1 | tr -d ' ')"
                ;;
        esac
    fi
done < "$CONFIG_FILE"

# Build RAW_PATH_KEYS/VALS: raw_sample name -> input path (before dedup)
declare -a RAW_PATH_KEYS=()
declare -a RAW_PATH_VALS=()
for i in "${!TREATMENTS[@]}"; do
    RAW_PATH_KEYS+=("${TREATMENTS[$i]}")
    RAW_PATH_VALS+=("${TREATMENT_PATHS[$i]}")
done
for i in "${!CONTROLS[@]}"; do
    RAW_PATH_KEYS+=("${CONTROLS[$i]}")
    RAW_PATH_VALS+=("${CONTROL_PATHS[$i]}")
done

# Collision check: same raw_sample name from two different paths is an error
for i in "${!RAW_PATH_KEYS[@]}"; do
    for j in "${!RAW_PATH_KEYS[@]}"; do
        if [[ "$i" -lt "$j" && "${RAW_PATH_KEYS[$i]}" == "${RAW_PATH_KEYS[$j]}" ]]; then
            if [[ "${RAW_PATH_VALS[$i]}" != "${RAW_PATH_VALS[$j]}" ]]; then
                echo "Error: raw_sample name '${RAW_PATH_KEYS[$i]}' is derived from two different paths:" >&2
                echo "  ${RAW_PATH_VALS[$i]}" >&2
                echo "  ${RAW_PATH_VALS[$j]}" >&2
                echo "Rename one of the input files or directories to resolve the conflict." >&2
                exit 1
            fi
        fi
    done
done

# Deduplicate RAW_PATH_KEYS/VALS (keep first occurrence of each name)
declare -a UNIQUE_PATH_KEYS=()
declare -a UNIQUE_PATH_VALS=()
for i in "${!RAW_PATH_KEYS[@]}"; do
    found="false"
    for k in "${UNIQUE_PATH_KEYS[@]+"${UNIQUE_PATH_KEYS[@]}"}"; do
        if [[ "$k" == "${RAW_PATH_KEYS[$i]}" ]]; then
            found="true"
            break
        fi
    done
    if [[ "$found" == "false" ]]; then
        UNIQUE_PATH_KEYS+=("${RAW_PATH_KEYS[$i]}")
        UNIQUE_PATH_VALS+=("${RAW_PATH_VALS[$i]}")
    fi
done

# Collect all raw samples (unique treatment + control names)
declare -a RAW_SAMPLES=()
for k in "${UNIQUE_PATH_KEYS[@]+"${UNIQUE_PATH_KEYS[@]}"}"; do
    RAW_SAMPLES+=("$k")
done
# Keep sorted for deterministic Snakefile output
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

if [[ ${#FEATURES[@]} -gt 0 ]]; then
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
fi

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
    if [[ $# -eq 0 ]]; then
        echo "[]"
        return
    fi
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
    if [[ ${#FEATURES[@]} -gt 0 ]]; then
        echo "FEATURES = $(python_list "${FEATURES[@]}")"
    else
        echo "FEATURES = []"
    fi
    echo ""
    echo "# Feature annotation parameters: feature -> (n_windows, window_size)"
    echo "FEATURE_PARAMS = {"
    if [[ ${#FEATURES[@]} -gt 0 ]]; then
        for i in "${!FEATURES[@]}"; do
            echo "    \"${FEATURES[$i]}\": (${FEATURE_N_WINDOWS[$i]}, ${FEATURE_WINDOW_SIZES[$i]}),"
        done
    fi
    echo "}"
    echo ""
    echo "# Window sizes for density calculation"
    if [[ ${#WINDOW_SIZES[@]} -gt 0 ]]; then
        echo "WINDOW_SIZES = $(python_list "${WINDOW_SIZES[@]}")"
    else
        echo "WINDOW_SIZES = []"
    fi
    echo ""
    echo "# Window sizes for mean reactivity"
    if [[ ${#MEAN_WINDOW_SIZES[@]} -gt 0 ]]; then
        echo "MEAN_WINDOW_SIZES = $(python_list "${MEAN_WINDOW_SIZES[@]}")"
    else
        echo "MEAN_WINDOW_SIZES = []"
    fi
    echo ""
    echo "# Significance threshold levels"
    if [[ ${#SIG_LEVELS[@]} -gt 0 ]]; then
        echo "SIG_LEVELS = $(python_list "${SIG_LEVELS[@]}")"
    else
        echo "SIG_LEVELS = []"
    fi
    echo ""
    echo "# Significance levels for density/bigwig filtering (sig=0 excluded from density)"
    echo "DENSITY_SIG_LEVELS = [s for s in SIG_LEVELS if s != \"0\"]"
    echo ""
    echo "# Strand directions"
    echo "STRANDS = [\"for\", \"rev\"]"
    echo ""
    echo "# All raw samples (treatment and control)"
    echo "RAW_SAMPLES = $(python_list "${RAW_SAMPLES[@]}")"
    echo ""
    echo "# Mapping from raw_sample name to input file or directory path"
    echo "RAW_PATHS = {"
    for i in "${!UNIQUE_PATH_KEYS[@]}"; do
        echo "    \"${UNIQUE_PATH_KEYS[$i]}\": \"${UNIQUE_PATH_VALS[$i]}\","
    done
    echo "}"
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
    if [[ ${#FEATURE_CHR_NAMES[@]} -gt 0 ]]; then
        for i in "${!FEATURE_CHR_NAMES[@]}"; do
            feature="${FEATURE_CHR_NAMES[$i]}"
            chrs_array=(${FEATURE_CHR_VALUES[$i]})
            echo "    \"$feature\": $(python_list "${chrs_array[@]}"),"
        done
    fi
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
    echo "    \"bg_mean\": \"data/bg_mean\" in TEMP_DIRS,"
    echo "    \"signal_reactivity\": \"data/signal_reactivity\" in TEMP_DIRS,"
    echo "    \"signal_bg\": \"data/signal_bg\" in TEMP_DIRS,"
    echo "    \"signal_bw\": \"data/signal_bw\" in TEMP_DIRS,"
    echo "    \"signal_bg_mean\": \"data/signal_bg_mean\" in TEMP_DIRS,"
    echo "}"
    echo ""
    echo "# Dorado basecalling model (used when pod5 input is detected)"
    echo "DORADO_MODEL = \"${DORADO_MODEL}\""
    echo ""
    echo "# IGV export settings"
    if [[ "$IGV_BAM" == "true" ]]; then
        echo "IGV_BAM = True"
    else
        echo "IGV_BAM = False"
    fi
    if [[ "$IGV_BIGWIG" == "true" ]]; then
        echo "IGV_BIGWIG = True"
    else
        echo "IGV_BIGWIG = False"
    fi
    echo "IGV_SOURCES = [\"all_alignments\", \"filtered_alignments\"]"
    echo ""
    echo "# Mapping from IGV source to input directory"
    echo "IGV_SOURCE_DIRS = {"
    echo "    \"all_alignments\": \"data/aligned_reads\","
    echo "    \"filtered_alignments\": \"data/filtered_alignments\","
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

include: "workflow/rules/phase0_pod5_processing.smk"
include: "workflow/rules/phase1_mapping_qc.smk"
include: "workflow/rules/phase2_perbase_error.smk"
include: "workflow/rules/phase2b_signal_deviation.smk"
include: "workflow/rules/phase3_reactivity.smk"
include: "workflow/rules/phase3b_signal_reactivity.smk"
include: "workflow/rules/phase4_analysis.smk"
include: "workflow/rules/phase4b_signal_analysis.smk"
include: "workflow/rules/phase5_annotate_features.smk"
include: "workflow/rules/phase6_igv.smk"
include: "workflow/rules/phase7_summary_tables_plots.smk"
include: "genome_specific_rules.smk"

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

        # Phase 3: Reactivity files (per chromosome)
        [f"data/reactivity/{genome}/{sample}_{strand}_{chr}.txt.gz"
         for genome in GENOMES
         for sample in SAMPLES
         for strand in STRANDS
         for chr in CHROMOSOMES[genome]],

        # Phase 4: Merged bigWig files
        expand("data/bw_merged/{genome}/significance_threshold_{sig}/{sample}_{strand}.bw",
               genome=GENOMES, sig=SIG_LEVELS, sample=SAMPLES, strand=STRANDS),

        # Phase 4: Merged density files (sig=0 excluded; density requires a significance threshold)
        expand("data/windows_merged/{genome}/window_size_{size}/significance_threshold_{sig}/{sample}_{strand}.bg",
               genome=GENOMES, size=WINDOW_SIZES, sig=DENSITY_SIG_LEVELS, sample=SAMPLES, strand=STRANDS)
        if WINDOW_SIZES and DENSITY_SIG_LEVELS else [],

        # Phase 4: Mean reactivity bigWig files
        expand("data/bw_mean_merged/{genome}/window_size_{mean_size}/{sample}_{strand}.bw",
               genome=GENOMES, mean_size=MEAN_WINDOW_SIZES, sample=SAMPLES, strand=STRANDS)
        if MEAN_WINDOW_SIZES else [],

        # Phase 5: Averaged feature annotations
        expand("data/annotations_averaged/{genome}/{feature}/{sample}_{strand}.txt.gz",
               genome=GENOMES, feature=FEATURES, sample=SAMPLES, strand=STRANDS),

        # Phase 6: IGV strand-split BAMs + indices
        expand("results/igv/{igv_source}/{genome}/{raw_sample}_{strand}.bam",
               igv_source=IGV_SOURCES, genome=GENOMES, raw_sample=RAW_SAMPLES, strand=STRANDS) +
        expand("results/igv/{igv_source}/{genome}/{raw_sample}_{strand}.bam.bai",
               igv_source=IGV_SOURCES, genome=GENOMES, raw_sample=RAW_SAMPLES, strand=STRANDS)
        if IGV_BAM else [],

        # Phase 6: IGV coverage bigWig files
        expand("results/igv/{igv_source}/{genome}/{raw_sample}_{strand}.bw",
               igv_source=IGV_SOURCES, genome=GENOMES, raw_sample=RAW_SAMPLES, strand=STRANDS)
        if IGV_BIGWIG else [],

        # Phase 7: Read statistics table
        "tables/read_stats_table.csv",

        # Phase 7: Alignment statistics table
        "tables/alignment_stats_table.csv",

        # Phase 7: Summary histograms (raw and filtered alignments, per genome per sample)
        expand("plots/histograms/{alignment}/{genome}/{raw_sample}_histograms.pdf",
               alignment=["aligned_reads", "filtered_alignments"],
               genome=GENOMES, raw_sample=RAW_SAMPLES),

        # Phase 7: Pairwise per-base error correlation tables and heatmaps (per genome)
        expand("tables/perbase_error_correlation/{genome}/Pairwise_correlation_table.csv",
               genome=GENOMES),
        expand("plots/perbase_error_correlation/{genome}/Pairwise_correlation_heatmap.pdf",
               genome=GENOMES),

        # Phase 7: Feature annotation plots (per genome per feature per sample)
        expand("plots/annotations_averaged/{genome}/{feature}/{sample}.pdf",
               genome=GENOMES, feature=FEATURES, sample=SAMPLES)
        if FEATURES else [],

        # Phase 2b: Per-base signal deviation (only for samples with pod5 input)
        [f"data/perbase_signal/{genome}/{rs}_{strand}.txt.gz"
         for genome in GENOMES
         for rs in RAW_SAMPLES
         for strand in STRANDS
         if has_pod5(rs)],

        # Phase 3b: Signal reactivity bigWig files (only for samples with pod5 input)
        [f"data/signal_bw_merged/{genome}/significance_threshold_{sig}/{sample}_{strand}.bw"
         for genome in GENOMES
         for sig in SIG_LEVELS
         for sample in SAMPLES
         for strand in STRANDS
         if has_pod5(get_treatment(sample)) and has_pod5(get_control(sample))],

        # Phase 4b: Mean signal reactivity bigWig (only when pod5 + mean windows configured)
        [f"data/signal_bw_mean_merged/{genome}/window_size_{mean_size}/{sample}_{strand}.bw"
         for genome in GENOMES
         for mean_size in MEAN_WINDOW_SIZES
         for sample in SAMPLES
         for strand in STRANDS
         if MEAN_WINDOW_SIZES and has_pod5(get_treatment(sample)) and has_pod5(get_control(sample))],

EOF

# ============================================================================
# Generate genome_specific_rules.smk
# ============================================================================
GENOME_RULES_FILE="genome_specific_rules.smk"

cat > "$GENOME_RULES_FILE" << 'EOF'
# ============================================================================
# Auto-generated genome-specific rules by CONFIG.sh
# Do not edit manually - regenerate from CONFIG file
# ============================================================================

EOF

# Generate split_perbase_by_chr and split_signal_by_chr rules for each genome
for i in "${!GENOME_CHR_NAMES[@]}"; do
    genome="${GENOME_CHR_NAMES[$i]}"
    chrs_array=(${GENOME_CHR_VALUES[$i]})
    chrs_python=$(python_list "${chrs_array[@]}")

    cat >> "$GENOME_RULES_FILE" << EOF
rule split_perbase_by_chr_${genome//./_}:
    """Split per-base error file by chromosome for ${genome}."""
    input:
        error="data/perbase_error/${genome}/{raw_sample}_{strand}.txt.gz"
    output:
        [wrap_output("perbase_error_by_chr", f) for f in
         expand("data/perbase_error_by_chr/${genome}/{{raw_sample}}_{{strand}}/{{raw_sample}}_{{strand}}_{chr}.txt.gz",
                chr=${chrs_python})]
    params:
        outdir="data/perbase_error_by_chr/${genome}/{raw_sample}_{strand}"
    log:
        "logs/split_perbase_by_chr/${genome}/{raw_sample}_{strand}.log"
    benchmark:
        "benchmarks/phase2/split_perbase_by_chr/${genome}/{raw_sample}_{strand}.tsv"
    wildcard_constraints:
        strand="for|rev"
    shell:
        """
        python3 workflow/scripts/Split_by_chr.sh \\
            -i {input.error} \\
            -d {params.outdir} \\
            2>&1 | tee {log}

        # Create empty gzipped files for any expected chromosomes with no data
        for f in {output}; do
            if [[ ! -f "\$f" ]]; then
                echo "No data for \$f — creating empty file" | tee -a {log}
                echo -n | gzip > "\$f"
            fi
        done

        echo "Successfully created chromosome files:" | tee -a {log}
        ls -la {params.outdir}/*.txt.gz | tee -a {log}
        """


rule split_signal_by_chr_${genome//./_}:
    """Split per-base signal deviation file by chromosome for ${genome}."""
    input:
        dev="data/perbase_signal/${genome}/{raw_sample}_{strand}.txt.gz"
    output:
        expand("data/perbase_signal_by_chr/${genome}/{{raw_sample}}_{{strand}}/{{raw_sample}}_{{strand}}_{chr}.txt.gz",
               chr=${chrs_python})
    params:
        outdir="data/perbase_signal_by_chr/${genome}/{raw_sample}_{strand}"
    log:
        "logs/split_signal_by_chr/${genome}/{raw_sample}_{strand}.log"
    benchmark:
        "benchmarks/phase2b/split_signal_by_chr/${genome}/{raw_sample}_{strand}.tsv"
    wildcard_constraints:
        strand="for|rev"
    shell:
        """
        python3 workflow/scripts/Split_by_chr.sh \\
            -i {input.dev} \\
            -d {params.outdir} \\
            2>&1 | tee {log}

        # Create empty gzipped files for any expected chromosomes with no data
        for f in {output}; do
            if [[ ! -f "\$f" ]]; then
                echo "No data for \$f — creating empty file" | tee -a {log}
                echo -n | gzip > "\$f"
            fi
        done

        echo "Successfully created chromosome files:" | tee -a {log}
        ls -la {params.outdir}/*.txt.gz | tee -a {log}
        """


EOF
done

echo "Generated $GENOME_RULES_FILE"
echo "Generated $OUTPUT_FILE from $CONFIG_FILE"
echo ""
echo "Summary:"
echo "  Genomes:        ${#GENOMES[@]} (${GENOMES[*]})"
echo "  Features:       ${#FEATURES[@]}$(if [[ ${#FEATURES[@]} -gt 0 ]]; then echo " (${FEATURES[*]})"; fi)"
echo "  Window sizes:   ${#WINDOW_SIZES[@]}$(if [[ ${#WINDOW_SIZES[@]} -gt 0 ]]; then echo " (${WINDOW_SIZES[*]})"; fi)"
echo "  Sig levels:     ${#SIG_LEVELS[@]}$(if [[ ${#SIG_LEVELS[@]} -gt 0 ]]; then echo " (${SIG_LEVELS[*]})"; fi)"
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
if [[ ${#FEATURE_CHR_NAMES[@]} -gt 0 ]]; then
    for i in "${!FEATURE_CHR_NAMES[@]}"; do
        feature="${FEATURE_CHR_NAMES[$i]}"
        chrs_array=(${FEATURE_CHR_VALUES[$i]})
        echo "  $feature: ${#chrs_array[@]} chromosomes"
    done
else
    echo "  (none)"
fi
echo ""
echo "Temporary directories:"
if [[ ${#TEMP_DIRS[@]} -gt 0 ]]; then
    for dir in "${TEMP_DIRS[@]}"; do
        echo "  $dir"
    done
else
    echo "  (none)"
fi
echo ""
echo "Dorado model (for pod5 inputs): $DORADO_MODEL"
echo "Raw input paths:"
for i in "${!UNIQUE_PATH_KEYS[@]}"; do
    echo "  ${UNIQUE_PATH_KEYS[$i]}: ${UNIQUE_PATH_VALS[$i]}"
done
