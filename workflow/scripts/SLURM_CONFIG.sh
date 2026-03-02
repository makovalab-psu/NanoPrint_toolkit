#!/bin/bash

# SLURM_CONFIG.sh - Generate a SLURM batch script for the NanoPrint pipeline
#
# Usage: SLURM_CONFIG.sh --env <conda_env> [-i Snakefile] [-o output.sh] [--alloc open] [--cores 8] [--slog dir]
#
# Runs a snakemake dry-run to count jobs, estimates wall time from benchmark
# data, and generates a ready-to-submit sbatch script for the Roar cluster.

set -euo pipefail

# ============================================================================
# Defaults
# ============================================================================
DATE_STAMP=$(date +%Y%m%d)
SNAKEFILE="./Snakefile"
OUTPUT_SCRIPT="${DATE_STAMP}_submit_nanoprint.sh"
ALLOC="open"
CORES=8
SLOG_DIR="${DATE_STAMP}_slurm_logs"
CONDA_ENV=""

# Locate benchmark summary relative to this script
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BENCHMARK_FILE="${SCRIPT_DIR}/../run_data/benchmark_summary.txt"

# ============================================================================
# Usage
# ============================================================================
usage() {
    cat <<EOF
Usage: $(basename "$0") --env <conda_env> [OPTIONS]

Generate a SLURM batch script for the NanoPrint pipeline on the Roar cluster.

Required:
  --env <path>       Path or name of the conda environment

Options:
  -i <file>          Input Snakefile (default: ./Snakefile)
  -o <file>          Output sbatch script (default: YYYYMMDD_submit_nanoprint.sh)
  --alloc <id>       Allocation (default: open)
                       open  -> --partition=open (free queue)
                       <id>  -> --partition=sla-prio --account=<id>
  --cores <n>        Number of cores (default: 8)
                       Each core gets 8 GB memory on the standard partition
  --slog <dir>       Directory for SLURM logs (default: YYYYMMDD_slurm_logs)
  -h, --help         Show this help message

Workflow:
  1. Runs snakemake --dry-run to count pending jobs
  2. Estimates wall time using benchmark data from previous runs
  3. Generates an sbatch script ready for submission

Example:
  $(basename "$0") --env nanoprint -i Snakefile --cores 16 --alloc open
  sbatch ${DATE_STAMP}_submit_nanoprint.sh
EOF
    exit 0
}

# ============================================================================
# Parse arguments (manual loop for long option support + bash 3.2 compat)
# ============================================================================
while [[ $# -gt 0 ]]; do
    case "$1" in
        --env)
            CONDA_ENV="$2"; shift 2 ;;
        -i)
            SNAKEFILE="$2"; shift 2 ;;
        -o)
            OUTPUT_SCRIPT="$2"; shift 2 ;;
        --alloc)
            ALLOC="$2"; shift 2 ;;
        --cores)
            CORES="$2"; shift 2 ;;
        --slog)
            SLOG_DIR="$2"; shift 2 ;;
        -h|--help)
            usage ;;
        *)
            echo "Error: Unknown option: $1" >&2
            echo "Run $(basename "$0") --help for usage." >&2
            exit 1 ;;
    esac
done

# ============================================================================
# Validate inputs
# ============================================================================
if [[ -z "$CONDA_ENV" ]]; then
    echo "Error: --env is required. Specify the conda environment name or path." >&2
    echo "Run $(basename "$0") --help for usage." >&2
    exit 1
fi

if [[ ! -f "$SNAKEFILE" ]]; then
    echo "Error: Snakefile not found: $SNAKEFILE" >&2
    exit 1
fi

if [[ ! -f "$BENCHMARK_FILE" ]]; then
    echo "Warning: Benchmark file not found: $BENCHMARK_FILE" >&2
    echo "         Time estimation will be skipped. Using default of 48 hours." >&2
    BENCHMARK_FILE=""
fi

# ============================================================================
# Run snakemake dry-run to get job summary
# ============================================================================
echo "Running snakemake dry-run..."

DRY_RUN_OUTPUT=$(snakemake --snakefile "$SNAKEFILE" --dry-run --quiet 2>&1) || {
    echo "Error: snakemake dry-run failed:" >&2
    echo "$DRY_RUN_OUTPUT" >&2
    exit 1
}

# Parse the job count summary table from dry-run output.
# Snakemake --dry-run --quiet prints a table like:
#   Job stats:
#   rule                      count
#   -------------------------  -----
#   alignment_stats               2
#   annotate_features             8
#   ...
#   total                        95
#
# We extract rule name and count, skipping header/separator/total lines.

# Store rule names and counts in parallel arrays (bash 3.2 compatible)
declare -a DRY_RULES=()
declare -a DRY_COUNTS=()
TOTAL_JOBS=0

while IFS= read -r line; do
    # Skip empty, header, separator, and total lines
    [[ -z "$line" ]] && continue
    [[ "$line" =~ ^[[:space:]]*Job ]] && continue
    [[ "$line" =~ ^[[:space:]]*rule ]] && continue
    [[ "$line" =~ ^[[:space:]]*--- ]] && continue
    [[ "$line" =~ ^[[:space:]]*total ]] && continue

    # Extract rule name and count
    rule=$(echo "$line" | awk '{print $1}')
    count=$(echo "$line" | awk '{print $NF}')

    # Validate we got a number
    if [[ "$count" =~ ^[0-9]+$ ]]; then
        DRY_RULES+=("$rule")
        DRY_COUNTS+=("$count")
        TOTAL_JOBS=$((TOTAL_JOBS + count))
    fi
done <<< "$DRY_RUN_OUTPUT"

echo "  Found $TOTAL_JOBS jobs across ${#DRY_RULES[@]} rules"

# ============================================================================
# Parse benchmark data and estimate wall time
# ============================================================================

# Convert human-readable time strings (e.g., "39.1m", "3.0h", "15.3s") to seconds
time_to_seconds() {
    local val="$1"
    if [[ "$val" == "-" ]] || [[ -z "$val" ]]; then
        echo "0"
        return
    fi
    # Strip any whitespace
    val=$(echo "$val" | tr -d ' ')
    local number="${val%[hms]}"
    local unit="${val: -1}"
    case "$unit" in
        h) echo "$number * 3600" | bc ;;
        m) echo "$number * 60" | bc ;;
        s) echo "$number" | bc ;;
        *)
            # Try as raw seconds
            echo "$number" | bc 2>/dev/null || echo "0"
            ;;
    esac
}

# Convert RSS strings (e.g., "6953MB") to integer MB
rss_to_mb() {
    local val="$1"
    if [[ "$val" == "-" ]] || [[ -z "$val" ]]; then
        echo "0"
        return
    fi
    echo "${val%MB}" | tr -d ' ' | bc 2>/dev/null || echo "0"
}

ESTIMATED_SECONDS=0
MAX_RSS_MB=0
MAX_RSS_RULE=""

if [[ -n "$BENCHMARK_FILE" ]] && [[ ${#DRY_RULES[@]} -gt 0 ]]; then
    echo ""
    echo "Estimating wall time from benchmark data..."
    echo ""
    printf "  %-28s %6s  %12s  %12s  %10s\n" "Rule" "Jobs" "Avg/job" "Subtotal" "PeakRSS"
    printf "  %-28s %6s  %12s  %12s  %10s\n" "----------------------------" "------" "------------" "------------" "----------"

    # For each rule in the dry-run, look up benchmark averages
    for i in "${!DRY_RULES[@]}"; do
        rule="${DRY_RULES[$i]}"
        count="${DRY_COUNTS[$i]}"

        # Extract bench03 and bench04 columns for this rule from the benchmark file.
        # Table columns (pipe-delimited, 1-indexed with leading empty field):
        #   2=Rule  3=bench03 Total  4=n  5=bench03 Avg  6=bench03 PeakRSS  7=bench03 AvgLoad
        #           8=bench04 Total  9=n 10=bench04 Avg 11=bench04 PeakRSS 12=bench04 AvgLoad
        bench_line=$(grep "^| ${rule} " "$BENCHMARK_FILE" 2>/dev/null || grep "^| ${rule}[[:space:]]" "$BENCHMARK_FILE" 2>/dev/null || echo "")

        avg_seconds=0
        avg_display="unknown"
        peak_rss_display="-"

        if [[ -n "$bench_line" ]]; then
            # Parse pipe-delimited columns
            bench03_avg=$(echo "$bench_line" | awk -F'|' '{print $5}' | tr -d ' ')
            bench03_rss=$(echo "$bench_line" | awk -F'|' '{print $6}' | tr -d ' ')
            bench04_avg=$(echo "$bench_line" | awk -F'|' '{print $10}' | tr -d ' ')
            bench04_rss=$(echo "$bench_line" | awk -F'|' '{print $11}' | tr -d ' ')

            bench03_sec=$(time_to_seconds "$bench03_avg")
            bench04_sec=$(time_to_seconds "$bench04_avg")

            # Use the maximum of the two as conservative estimate
            if (( $(echo "$bench03_sec > $bench04_sec" | bc -l) )); then
                avg_seconds="$bench03_sec"
                avg_display="$bench03_avg"
            else
                avg_seconds="$bench04_sec"
                avg_display="$bench04_avg"
            fi

            # Track peak RSS: use the larger of bench03/bench04
            bench03_rss_mb=$(rss_to_mb "$bench03_rss")
            bench04_rss_mb=$(rss_to_mb "$bench04_rss")
            if (( $(echo "$bench03_rss_mb > $bench04_rss_mb" | bc -l) )); then
                rule_rss_mb="$bench03_rss_mb"
                peak_rss_display="$bench03_rss"
            else
                rule_rss_mb="$bench04_rss_mb"
                peak_rss_display="$bench04_rss"
            fi

            if (( $(echo "$rule_rss_mb > $MAX_RSS_MB" | bc -l) )); then
                MAX_RSS_MB="$rule_rss_mb"
                MAX_RSS_RULE="$rule"
            fi
        fi

        subtotal=$(echo "$avg_seconds * $count" | bc)
        ESTIMATED_SECONDS=$(echo "$ESTIMATED_SECONDS + $subtotal" | bc)

        # Format subtotal for display
        sub_hours=$(echo "$subtotal / 3600" | bc)
        sub_remainder=$(echo "$subtotal - ($sub_hours * 3600)" | bc)
        sub_min=$(echo "$sub_remainder / 60" | bc)
        if [[ "$sub_hours" -gt 0 ]]; then
            sub_display="${sub_hours}h ${sub_min}m"
        elif [[ $(echo "$subtotal > 60" | bc) -eq 1 ]]; then
            sub_display="${sub_min}m"
        else
            sub_display="${subtotal}s"
        fi

        printf "  %-28s %6d  %12s  %12s  %10s\n" "$rule" "$count" "$avg_display" "$sub_display" "$peak_rss_display"
    done

    echo ""

    # Account for parallelism: divide by number of cores
    PARALLEL_SECONDS=$(echo "$ESTIMATED_SECONDS / $CORES" | bc)

    # Add 20% safety margin
    SAFE_SECONDS=$(echo "$PARALLEL_SECONDS * 1.2 / 1" | bc)

    # Round up to the next hour
    WALL_HOURS=$(( (SAFE_SECONDS + 3599) / 3600 ))

    # Cap at 14 days (336 hours) - Roar max for normal QOS
    if [[ $WALL_HOURS -gt 336 ]]; then
        echo "  Warning: Estimated time exceeds 14-day limit. Capping at 336 hours."
        WALL_HOURS=336
    fi

    # Minimum 1 hour
    if [[ $WALL_HOURS -lt 1 ]]; then
        WALL_HOURS=1
    fi

    serial_h=$(echo "$ESTIMATED_SECONDS / 3600" | bc)
    parallel_h=$(echo "$PARALLEL_SECONDS / 3600" | bc)
    echo "  Total serial time:    ~${serial_h} hours"
    echo "  Parallel (${CORES} cores): ~${parallel_h} hours"
    echo "  With 20% margin:     ${WALL_HOURS} hours"

    # Memory summary
    if [[ "$MAX_RSS_MB" -gt 0 ]]; then
        echo ""
        echo "Peak memory requirements:"
        # Compute minimum cores needed at 8 GB/core (8192 MB/core)
        min_cores_mem=$(echo "($MAX_RSS_MB + 8191) / 8192" | bc)
        [[ "$min_cores_mem" -lt 1 ]] && min_cores_mem=1
        echo "  Rule with max RSS:  ${MAX_RSS_RULE} (${MAX_RSS_MB}MB)"
        echo "  Min cores for mem:  ${min_cores_mem} (at 8GB/core on standard partition)"
        total_mem_gb=$(( CORES * 8 ))
        if [[ "$CORES" -lt "$min_cores_mem" ]]; then
            echo "  WARNING: --cores ${CORES} provides ${total_mem_gb}GB; ${MAX_RSS_RULE} may need $((min_cores_mem * 8))GB"
        else
            echo "  --cores ${CORES} provides ${total_mem_gb}GB -- sufficient"
        fi
    fi
else
    # No benchmarks available - use a safe default
    WALL_HOURS=48
    echo "  Using default wall time: ${WALL_HOURS} hours"
fi

# Format time as HH:MM:SS for SBATCH
WALL_TIME="${WALL_HOURS}:00:00"

# ============================================================================
# Determine partition and account settings
# ============================================================================
if [[ "$ALLOC" == "open" ]]; then
    PARTITION_LINE="#SBATCH --partition=open"
    ACCOUNT_LINE=""
else
    PARTITION_LINE="#SBATCH --partition=sla-prio"
    ACCOUNT_LINE="#SBATCH --account=${ALLOC}"
fi

# Memory: standard partition = 8 GB/core
MEM_PER_CPU="8gb"

# ============================================================================
# Build SBATCH header
# ============================================================================
SBATCH_HEADER="#SBATCH --job-name=nanoprint
${PARTITION_LINE}"
if [[ -n "$ACCOUNT_LINE" ]]; then
    SBATCH_HEADER="${SBATCH_HEADER}
${ACCOUNT_LINE}"
fi
SBATCH_HEADER="${SBATCH_HEADER}
#SBATCH --nodes=1
#SBATCH --ntasks=${CORES}
#SBATCH --mem-per-cpu=${MEM_PER_CPU}
#SBATCH --time=${WALL_TIME}
#SBATCH --output=${SLOG_DIR}/nanoprint_%j.out
#SBATCH --error=${SLOG_DIR}/nanoprint_%j.err"

# ============================================================================
# Generate the sbatch script
# ============================================================================
echo ""
echo "Generating sbatch script: $OUTPUT_SCRIPT"

cat > "$OUTPUT_SCRIPT" << SBATCH_EOF
#!/bin/bash
# ============================================================================
# NanoPrint Pipeline - SLURM Batch Script
# Generated by SLURM_CONFIG.sh on $(date '+%Y-%m-%d %H:%M:%S')
# ============================================================================

${SBATCH_HEADER}

# ============================================================================
# Setup
# ============================================================================

# Create log directory if it does not exist
mkdir -p ${SLOG_DIR}

# Print job info
echo "============================================"
echo "NanoPrint Pipeline - SLURM Job"
echo "============================================"
echo "Job ID:       \$SLURM_JOB_ID"
echo "Job Name:     \$SLURM_JOB_NAME"
echo "Nodes:        \$SLURM_NODELIST"
echo "Cores:        \$SLURM_NTASKS"
echo "Partition:    \$SLURM_QUEUE"
echo "Submit Dir:   \$SLURM_SUBMIT_DIR"
echo "Start Time:   \$(date)"
echo "============================================"
echo ""

# ============================================================================
# Load software
# ============================================================================
module load anaconda
conda activate ${CONDA_ENV}

# ============================================================================
# Run pipeline
# ============================================================================
echo "Starting snakemake pipeline..."
echo ""

snakemake \\
    --snakefile ${SNAKEFILE} \\
    --cores \$SLURM_NTASKS

echo ""
echo "============================================"
echo "Pipeline finished: \$(date)"
echo "============================================"

# ============================================================================
# Report resource usage
# ============================================================================
echo ""
echo "Resource usage:"
sacct -j \$SLURM_JOB_ID --format=JobID,JobName,MaxRSS,Elapsed,TotalCPU,State
SBATCH_EOF

chmod +x "$OUTPUT_SCRIPT"

# ============================================================================
# Print summary
# ============================================================================
echo ""
echo "============================================"
echo "  SLURM Configuration Summary"
echo "============================================"
echo "  Snakefile:    $SNAKEFILE"
echo "  Conda env:    $CONDA_ENV"
echo "  Allocation:   $ALLOC"
if [[ "$ALLOC" == "open" ]]; then
    echo "  Partition:    open"
else
    echo "  Partition:    sla-prio"
fi
echo "  Cores:        $CORES"
echo "  Memory:       ${MEM_PER_CPU}/core ($(( CORES * 8 )) GB total)"
echo "  Wall time:    $WALL_TIME"
echo "  Log dir:      $SLOG_DIR"
echo "  Total jobs:   $TOTAL_JOBS"
echo "  Output:       $OUTPUT_SCRIPT"
echo "============================================"
echo ""
echo "To submit:"
echo "  sbatch $OUTPUT_SCRIPT"
echo ""
