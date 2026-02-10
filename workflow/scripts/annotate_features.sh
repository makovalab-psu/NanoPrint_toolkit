#!/bin/bash

# Annotate Features: Calculate signal around genomic features
# Works on per-chromosome split files
# Uses chunk-based processing to limit memory usage

set -euo pipefail

# Default parameters
N_WINDOWS=1000
WINDOW_SIZE=10
TMP_DIR=""
CHUNK_SIZE=100000

# Usage function
usage() {
    cat << EOF
Usage: $(basename "$0") -b <annotations.bed> -t <treatment.txt> -c <control.txt> -r <reactivity.txt> -g <genome.fai> -o <output.txt.gz> [OPTIONS]

Calculate functional genomic signal surrounding genomic features.

Required arguments:
    -b    BED file with genomic features (split by chromosome)
    -t    Treatment per-base error file (split by chromosome)
    -c    Control per-base error file (split by chromosome)
    -r    Reactivity file (split by chromosome)
    -g    Genome index file (.fai) for chromosome sizes
    -o    Output file (gzipped)

Optional arguments:
    -n, --n-windows NUM     Number of windows on each side (default: 1000)
    -w, --window-size NUM   Window size in nucleotides (default: 10)
    -T DIRECTORY            Temporary directory (default: output directory)
    -h                      Show this help message

Input formats:
    BED file: Standard BED6 format (chr, start, end, name, score, strand)
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
    $(basename "$0") -b features_chr1.bed -t treat_for_chr1.txt -c ctrl_for_chr1.txt -r react_for_chr1.txt -g genome.fa.fai -o output_chr1.txt.gz
EOF
    exit 1
}

# Parse command line arguments
BED_FILE=""
TREATMENT=""
CONTROL=""
REACTIVITY=""
GENOME_FAI=""
OUTPUT=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -b) BED_FILE="$2"; shift 2 ;;
        -t) TREATMENT="$2"; shift 2 ;;
        -c) CONTROL="$2"; shift 2 ;;
        -r) REACTIVITY="$2"; shift 2 ;;
        -g) GENOME_FAI="$2"; shift 2 ;;
        -o) OUTPUT="$2"; shift 2 ;;
        -n|--n-windows) N_WINDOWS="$2"; shift 2 ;;
        -w|--window-size) WINDOW_SIZE="$2"; shift 2 ;;
        -T) TMP_DIR="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Error: Unknown option: $1" >&2; usage ;;
    esac
done

# Validate required arguments
if [[ -z "$BED_FILE" || -z "$TREATMENT" || -z "$CONTROL" || -z "$REACTIVITY" || -z "$GENOME_FAI" || -z "$OUTPUT" ]]; then
    echo "Error: Missing required arguments" >&2
    usage
fi

# Check input files exist
for f in "$BED_FILE" "$TREATMENT" "$CONTROL" "$REACTIVITY" "$GENOME_FAI"; do
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
echo "Genome FAI: $GENOME_FAI"
echo "Output: $OUTPUT"
echo "Windows: $N_WINDOWS on each side, ${WINDOW_SIZE}bp each"
echo "Chunk size: ${CHUNK_SIZE}bp"
echo ""

# Process features using Python with chunk-based memory management
echo "Processing features with chunk-based memory management..."
python3 - "$BED_FILE" "$TREATMENT" "$CONTROL" "$REACTIVITY" "$GENOME_FAI" "$OUTPUT" "$N_WINDOWS" "$WINDOW_SIZE" "$TMP_DIR" "$CHUNK_SIZE" << 'PYTHON_SCRIPT'
import sys
import os
import gzip

bed_file = sys.argv[1]
treatment_file = sys.argv[2]
control_file = sys.argv[3]
reactivity_file = sys.argv[4]
fai_file = sys.argv[5]
output_file = sys.argv[6]
n_windows = int(sys.argv[7])
window_size = int(sys.argv[8])
tmp_dir = sys.argv[9]
chunk_size = int(sys.argv[10])

window_span = n_windows * window_size

def open_file(filepath):
    """Open file, handling gzip if needed."""
    if filepath.endswith('.gz'):
        return gzip.open(filepath, 'rt')
    return open(filepath, 'r')

def get_chromosome_from_file(filepath):
    """Read first line to get chromosome name from column 1."""
    with open_file(filepath) as f:
        line = f.readline()
        if line:
            return line.strip().split('\t')[0]
    return None

def get_chrom_size_from_fai(fai_path, chrom):
    """Look up chromosome size from .fai index file."""
    with open(fai_path) as f:
        for line in f:
            fields = line.strip().split('\t')
            if fields[0] == chrom:
                return int(fields[1])
    raise ValueError(f"Chromosome {chrom} not found in {fai_path}")

def split_perbase_into_bins(filepath, out_prefix, chunk_size):
    """Split per-base file into bin files (single pass, no padding)."""
    os.makedirs(out_prefix, exist_ok=True)
    handles = {}
    with open_file(filepath) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue
            pos = int(fields[1])
            bin_idx = pos // chunk_size
            if bin_idx not in handles:
                handles[bin_idx] = open(f"{out_prefix}/bin_{bin_idx}.txt", 'w')
            handles[bin_idx].write(line)
    for h in handles.values():
        h.close()
    return set(handles.keys())

def split_reactivity_into_bins(filepath, out_prefix, chunk_size):
    """Split reactivity file into bin files (single pass, no padding)."""
    os.makedirs(out_prefix, exist_ok=True)
    handles = {}
    with open_file(filepath) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 4:
                continue
            pos = int(fields[1])
            bin_idx = pos // chunk_size
            if bin_idx not in handles:
                handles[bin_idx] = open(f"{out_prefix}/bin_{bin_idx}.txt", 'w')
            handles[bin_idx].write(line)
    for h in handles.values():
        h.close()
    return set(handles.keys())

def split_bed_into_bins(filepath, out_prefix, chunk_size):
    """Split BED file into bins by reference point."""
    os.makedirs(out_prefix, exist_ok=True)
    handles = {}
    with open_file(filepath) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            strand = fields[5] if len(fields) > 5 else '+'
            ref_pos = int(fields[2]) if strand == '-' else int(fields[1])
            bin_idx = ref_pos // chunk_size
            if bin_idx not in handles:
                handles[bin_idx] = open(f"{out_prefix}/bin_{bin_idx}.bed", 'w')
            handles[bin_idx].write(line)
    for h in handles.values():
        h.close()
    return set(handles.keys())

def load_perbase_bin(filepath):
    """Load per-base error bin file into dictionary keyed by position."""
    data = {}
    if not os.path.exists(filepath):
        return data
    with open(filepath, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 5:
                pos = int(fields[1])
                cov = float(fields[3])
                err = float(fields[4])
                data[pos] = {'coverage': cov, 'error': err}
    return data

def load_reactivity_bin(filepath):
    """Load reactivity bin file into dictionary keyed by position."""
    data = {}
    if not os.path.exists(filepath):
        return data
    with open(filepath, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 4:
                pos = int(fields[1])
                try:
                    react = float(fields[3])
                    if react == 999999 or react == -999999:
                        continue
                    data[pos] = react
                except ValueError:
                    continue
    return data

def load_bed_bin(filepath):
    """Load BED bin file into list of features."""
    features = []
    if not os.path.exists(filepath):
        return features
    with open(filepath, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            features.append({
                'start': int(fields[1]),
                'end': int(fields[2]),
                'strand': fields[5] if len(fields) > 5 else '+'
            })
    return features

class ChunkCache:
    """Cache that holds at most 2 adjacent bins loaded."""
    def __init__(self, bin_dir, loader_func, file_ext='txt'):
        self.bin_dir = bin_dir
        self.loader = loader_func
        self.file_ext = file_ext
        self.cache = {}  # bin_idx -> data dict

    def get_data(self, needed_bins):
        """Ensure needed bins are loaded, unload others."""
        needed = set(needed_bins)
        # Unload bins no longer needed
        for idx in list(self.cache.keys()):
            if idx not in needed:
                del self.cache[idx]
        # Load missing bins
        for idx in needed:
            if idx not in self.cache:
                filepath = f"{self.bin_dir}/bin_{idx}.{self.file_ext}"
                self.cache[idx] = self.loader(filepath)
        # Merge loaded bins into single dict for lookup
        merged = {}
        for data in self.cache.values():
            merged.update(data)
        return merged

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

# Get chromosome from first line of treatment file
print("Determining chromosome...", file=sys.stderr)
chrom = get_chromosome_from_file(treatment_file)
if chrom is None:
    print("Error: Could not determine chromosome from treatment file", file=sys.stderr)
    sys.exit(1)
print(f"  Chromosome: {chrom}", file=sys.stderr)

# Get chromosome size from FAI file
chrom_size = get_chrom_size_from_fai(fai_file, chrom)
print(f"  Chromosome size: {chrom_size:,}", file=sys.stderr)

num_bins = (chrom_size + chunk_size - 1) // chunk_size
print(f"  Number of bins: {num_bins}", file=sys.stderr)

# Determine strand from filename
strand = "for"
if "_rev_" in treatment_file or "_rev." in treatment_file:
    strand = "rev"

# Pass 1: Split all files into bin temp files
print("Pass 1: Splitting files into bins...", file=sys.stderr)

print("  Splitting treatment file...", file=sys.stderr)
treat_bins = split_perbase_into_bins(treatment_file, f"{tmp_dir}/treatment", chunk_size)
print(f"    Created {len(treat_bins)} bins", file=sys.stderr)

print("  Splitting control file...", file=sys.stderr)
ctrl_bins = split_perbase_into_bins(control_file, f"{tmp_dir}/control", chunk_size)
print(f"    Created {len(ctrl_bins)} bins", file=sys.stderr)

print("  Splitting reactivity file...", file=sys.stderr)
react_bins = split_reactivity_into_bins(reactivity_file, f"{tmp_dir}/reactivity", chunk_size)
print(f"    Created {len(react_bins)} bins", file=sys.stderr)

print("  Splitting BED file...", file=sys.stderr)
feature_bins = split_bed_into_bins(bed_file, f"{tmp_dir}/features", chunk_size)
print(f"    Created {len(feature_bins)} bins", file=sys.stderr)

# Pass 2: Process features with sliding window cache
print("Pass 2: Processing features with sliding cache...", file=sys.stderr)

treat_cache = ChunkCache(f"{tmp_dir}/treatment", load_perbase_bin)
ctrl_cache = ChunkCache(f"{tmp_dir}/control", load_perbase_bin)
react_cache = ChunkCache(f"{tmp_dir}/reactivity", load_reactivity_bin)

feature_count = 0
with gzip.open(output_file, 'wt') as out:
    # Write header
    out.write("Distance\tCoverage\tPerbase_error\tReactivity\tSample\tStrand\n")

    # Process each bin's features in order
    for bin_idx in sorted(feature_bins):
        features_file = f"{tmp_dir}/features/bin_{bin_idx}.bed"
        features = load_bed_bin(features_file)

        for feature in features:
            feature_strand = feature['strand']
            if feature_strand == '-':
                ref_pos = feature['end']
            else:
                ref_pos = feature['start']

            # Determine which bins we need for window span (max 2)
            window_start_pos = ref_pos - window_span
            window_end_pos = ref_pos + window_span
            needed_bins = set()
            needed_bins.add(max(0, window_start_pos // chunk_size))
            needed_bins.add(min(num_bins - 1, window_end_pos // chunk_size))

            # Load needed chunks (cache handles load/unload)
            treat_data = treat_cache.get_data(needed_bins)
            ctrl_data = ctrl_cache.get_data(needed_bins)
            react_data = react_cache.get_data(needed_bins)

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
                treat_cov = calculate_window_coverage(treat_data, window_start, window_end)
                treat_err = calculate_window_average(treat_data, window_start, window_end, 'error')
                react_val = calculate_window_average(react_data, window_start, window_end)

                if treat_cov is not None and treat_err is not None:
                    react_str = f"{react_val:.6f}" if react_val is not None else ""
                    out.write(f"{distance}\t{treat_cov:.2f}\t{treat_err:.6f}\t{react_str}\tTreatment\t{strand}\n")

                # Control
                ctrl_cov = calculate_window_coverage(ctrl_data, window_start, window_end)
                ctrl_err = calculate_window_average(ctrl_data, window_start, window_end, 'error')

                if ctrl_cov is not None and ctrl_err is not None:
                    out.write(f"{distance}\t{ctrl_cov:.2f}\t{ctrl_err:.6f}\t\tControl\t{strand}\n")

        if bin_idx % 100 == 0:
            print(f"  Processed bin {bin_idx}/{num_bins-1}...", file=sys.stderr)

print(f"Processed {feature_count} features", file=sys.stderr)
PYTHON_SCRIPT

echo ""
echo "Done. Output: $OUTPUT"
