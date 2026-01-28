#!/bin/bash

# Annotate Features: Calculate signal around genomic features
# Works on per-chromosome split files

set -euo pipefail

# Default parameters
N_WINDOWS=1000
WINDOW_SIZE=10
TMP_DIR=""

# Usage function
usage() {
    cat << EOF
Usage: $(basename "$0") -b <annotations.bed> -t <treatment.txt> -c <control.txt> -r <reactivity.txt> -o <output.txt.gz> [OPTIONS]

Calculate functional genomic signal surrounding genomic features.

Required arguments:
    -b    BED file with genomic features (split by chromosome)
    -t    Treatment per-base error file (split by chromosome)
    -c    Control per-base error file (split by chromosome)
    -r    Reactivity file (split by chromosome)
    -o    Output file (gzipped)

Optional arguments:
    -n, --n-windows NUM     Number of windows on each side (default: 1000)
    -w, --window-size NUM   Window size in nucleotides (default: 10)
    -T DIRECTORY            Temporary directory (default: output directory)
    -h                      Show this help message

Input formats:
    BED file: chr, start, end, strand (tab-delimited)
    Per-base error: chr, position, nucleotide, coverage, error (tab-delimited)
    Reactivity: chr, position, nucleotide, reactivity (tab-delimited)

Output format (tab-delimited, gzipped):
    Distance        - Distance from feature reference point
    Coverage        - Average coverage in window
    Perbase_error   - Average per-base error in window
    Reactivity      - Average reactivity in window (Treatment only)
    Sample          - Treatment or Control
    Strand          - for or rev

Example:
    $(basename "$0") -b features_chr1.bed -t treat_for_chr1.txt -c ctrl_for_chr1.txt -r react_for_chr1.txt -o output_chr1.txt.gz
EOF
    exit 1
}

# Parse command line arguments
BED_FILE=""
TREATMENT=""
CONTROL=""
REACTIVITY=""
OUTPUT=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -b) BED_FILE="$2"; shift 2 ;;
        -t) TREATMENT="$2"; shift 2 ;;
        -c) CONTROL="$2"; shift 2 ;;
        -r) REACTIVITY="$2"; shift 2 ;;
        -o) OUTPUT="$2"; shift 2 ;;
        -n|--n-windows) N_WINDOWS="$2"; shift 2 ;;
        -w|--window-size) WINDOW_SIZE="$2"; shift 2 ;;
        -T) TMP_DIR="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Error: Unknown option: $1" >&2; usage ;;
    esac
done

# Validate required arguments
if [[ -z "$BED_FILE" || -z "$TREATMENT" || -z "$CONTROL" || -z "$REACTIVITY" || -z "$OUTPUT" ]]; then
    echo "Error: Missing required arguments" >&2
    usage
fi

# Check input files exist
for f in "$BED_FILE" "$TREATMENT" "$CONTROL" "$REACTIVITY"; do
    if [[ ! -f "$f" ]]; then
        echo "Error: File not found: $f" >&2
        exit 1
    fi
done

# Get output directory
OUT_DIR=$(dirname "$OUTPUT")
if [[ -z "$OUT_DIR" || "$OUT_DIR" == "." ]]; then
    OUT_DIR="$(pwd)"
fi
mkdir -p "$OUT_DIR"

# Set temp directory
if [[ -z "$TMP_DIR" ]]; then
    TMP_DIR="${OUT_DIR}/tmp_annotate_$$"
fi
mkdir -p "$TMP_DIR"

# Cleanup function
cleanup() {
    if [[ -d "$TMP_DIR" ]]; then
        rm -rf "$TMP_DIR"
    fi
}
trap cleanup EXIT INT TERM

echo "=== Annotate Features ==="
echo "BED file: $BED_FILE"
echo "Treatment: $TREATMENT"
echo "Control: $CONTROL"
echo "Reactivity: $REACTIVITY"
echo "Output: $OUTPUT"
echo "Windows: $N_WINDOWS on each side, ${WINDOW_SIZE}bp each"
echo ""

# Process features using Python
echo "Processing features..."
python3 - "$BED_FILE" "$TREATMENT" "$CONTROL" "$REACTIVITY" "$OUTPUT" "$N_WINDOWS" "$WINDOW_SIZE" << 'PYTHON_SCRIPT'
import sys
import gzip
from collections import defaultdict

bed_file = sys.argv[1]
treatment_file = sys.argv[2]
control_file = sys.argv[3]
reactivity_file = sys.argv[4]
output_file = sys.argv[5]
n_windows = int(sys.argv[6])
window_size = int(sys.argv[7])

def load_perbase_file(filepath):
    """Load per-base error file into dictionary keyed by position."""
    data = {}
    with open(filepath, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 5:
                pos = int(fields[1])
                cov = float(fields[3])
                err = float(fields[4])
                data[pos] = {'coverage': cov, 'error': err}
    return data

def load_reactivity_file(filepath):
    """Load reactivity file into dictionary keyed by position."""
    data = {}
    with open(filepath, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 4:
                pos = int(fields[1])
                # Handle special codes
                try:
                    react = float(fields[3])
                    if react == 999999 or react == -999999:
                        continue
                    data[pos] = react
                except ValueError:
                    continue
    return data

def calculate_window_average(data_dict, start, end, key='error'):
    """Calculate average value in a window."""
    values = []
    for pos in range(start, end + 1):
        if pos in data_dict:
            if isinstance(data_dict[pos], dict):
                values.append(data_dict[pos][key])
            else:
                values.append(data_dict[pos])
    if values:
        return sum(values) / len(values)
    return None

def calculate_window_coverage(data_dict, start, end):
    """Calculate average coverage in a window."""
    coverages = []
    for pos in range(start, end + 1):
        if pos in data_dict and 'coverage' in data_dict[pos]:
            coverages.append(data_dict[pos]['coverage'])
    if coverages:
        return sum(coverages) / len(coverages)
    return None

print("Loading treatment data...", file=sys.stderr)
treatment_data = load_perbase_file(treatment_file)
print(f"  Loaded {len(treatment_data)} positions", file=sys.stderr)

print("Loading control data...", file=sys.stderr)
control_data = load_perbase_file(control_file)
print(f"  Loaded {len(control_data)} positions", file=sys.stderr)

print("Loading reactivity data...", file=sys.stderr)
reactivity_data = load_reactivity_file(reactivity_file)
print(f"  Loaded {len(reactivity_data)} positions", file=sys.stderr)

# Determine strand from filename (assumes format like *_for_*.txt or *_rev_*.txt)
strand = "for"
if "_rev_" in treatment_file or "_rev." in treatment_file:
    strand = "rev"

print(f"Processing features (strand: {strand})...", file=sys.stderr)

with gzip.open(output_file, 'wt') as out:
    # Write header
    out.write("Distance\tCoverage\tPerbase_error\tReactivity\tSample\tStrand\n")

    feature_count = 0
    with open(bed_file, 'r') as bed:
        for line in bed:
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue

            start = int(fields[1])
            end = int(fields[2])
            feature_strand = fields[3] if len(fields) > 3 else '+'

            # Determine reference position
            if feature_strand == '-':
                ref_pos = end
            else:
                ref_pos = start

            feature_count += 1

            # Process each window
            for i in range(-n_windows, n_windows + 1):
                # Calculate distance (flip for reverse strand features)
                if feature_strand == '-':
                    distance = -i * window_size
                else:
                    distance = i * window_size

                window_center = ref_pos + (i * window_size)
                window_start = window_center - window_size // 2
                window_end = window_center + window_size // 2

                # Treatment
                treat_cov = calculate_window_coverage(treatment_data, window_start, window_end)
                treat_err = calculate_window_average(treatment_data, window_start, window_end, 'error')
                react_val = calculate_window_average(reactivity_data, window_start, window_end)

                if treat_cov is not None and treat_err is not None:
                    react_str = f"{react_val:.6f}" if react_val is not None else ""
                    out.write(f"{distance}\t{treat_cov:.2f}\t{treat_err:.6f}\t{react_str}\tTreatment\t{strand}\n")

                # Control
                ctrl_cov = calculate_window_coverage(control_data, window_start, window_end)
                ctrl_err = calculate_window_average(control_data, window_start, window_end, 'error')

                if ctrl_cov is not None and ctrl_err is not None:
                    out.write(f"{distance}\t{ctrl_cov:.2f}\t{ctrl_err:.6f}\t\tControl\t{strand}\n")

print(f"Processed {feature_count} features", file=sys.stderr)
PYTHON_SCRIPT

echo ""
echo "Done. Output: $OUTPUT"
