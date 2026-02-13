# NanoPrint Toolkit - Workflow Summary

## Overview

NanoPrint_toolkit is a Snakemake-based pipeline for analyzing **chemical footprinting data** from Oxford Nanopore long-read sequencing. It processes raw sequencing reads through alignment, per-base error quantification, and reactivity calculation (treatment minus control).

## Key Concepts

- **Per-base error**: Error rate at each genomic position, calculated from aligned reads
- **Reactivity**: Difference between treatment and control per-base error rates
- **Parallelization**: Data is split by chromosome for efficient parallel processing
- **Raw samples**: Individual treatment/control FASTQ/BAM files processed in phases 1-2
- **Relationship samples**: Combined sample names used in phases 3-5 outputs

## Configuration

### CONFIG File Format

The pipeline is configured via a CONFIG file with prefix notation:

| Prefix | Description | Example |
|--------|-------------|---------|
| `^g` | Genome file (in resources/genomes/) | `^g test_genome.fa` |
| `^f` | Feature file (in resources/features/) | `^f g4Discovery.bed` |
| `^w` | Window size for density calculation | `^w 1000000` |
| `^s` | Significance threshold (1-4) | `^s 2` |
| `^r` | Relationship: Sample, Treatment, Control | `^r SampleName Treatment.bam Control.bam` |
| `^t` | Temporary directory (auto-deleted) | `^t data/aligned_reads` |

Extensions are automatically stripped from all file values.

### Generating the Snakefile

```bash
./workflow/scripts/CONFIG.sh -i CONFIG -o Snakefile
```

This parses CONFIG and generates a Snakefile with:
- Hardcoded wildcard lists (GENOMES, FEATURES, etc.)
- RELATIONSHIPS dict mapping sample → (treatment, control)
- `rule all` with expand() for all target outputs

## Wildcards

| Wildcard | Phase | Description |
|----------|-------|-------------|
| `{genome}` | All | Reference genome name |
| `{raw_sample}` | 1-2 | Individual treatment/control sample |
| `{raw_sample_a}`, `{raw_sample_b}` | 2 | Samples for correlation |
| `{sample}` | 3-5 | Relationship sample name (output prefix) |
| `{treatment_sample}`, `{control_sample}` | 3, 5 | Raw sample names for inputs |
| `{strand}` | 2-5 | `for` or `rev` |
| `{chr}` | 3-5 | Chromosome (from checkpoints) |
| `{size}` | 4 | Window size |
| `{sig}` | 4 | Significance threshold level |
| `{feature}` | 5 | Feature name |
| `{alignment}` | 1 | `aligned_reads` or `filtered_alignments` |

## Sample Flow

```
CONFIG Relationship:
^r  Hsap_HG002_LCL  Hsap_HG002_LCL_Mn04.bam  Hsap_HG002_LCL_CTRL.bam
    └─ {sample}     └─ {treatment_sample}     └─ {control_sample}
                    └────────────────────────────────────────────┘
                              These are {raw_sample} values

Phases 1-2: Process EACH raw sample independently
┌─────────────────────────────────────────────────────────────────┐
│  Hsap_HG002_LCL_Mn04 → perbase_error/Hsap_HG002_LCL_Mn04_for.txt.gz
│  Hsap_HG002_LCL_CTRL → perbase_error/Hsap_HG002_LCL_CTRL_for.txt.gz
└─────────────────────────────────────────────────────────────────┘

Phases 3-5: Combine using {sample} from relationship
┌─────────────────────────────────────────────────────────────────┐
│  treatment: Hsap_HG002_LCL_Mn04 ─┐
│                                  ├→ reactivity/Hsap_HG002_LCL_for_chr1.txt.gz
│  control:   Hsap_HG002_LCL_CTRL ─┘
└─────────────────────────────────────────────────────────────────┘
```

## Input File Flexibility

Phase 1 rules use `find_raw_reads()` helper function that automatically detects input format:
- Checks for `.fastq.gz`, `.fastq`, `.bam` in order
- No need to specify extension in CONFIG

## Pipeline Phases

### Phase 1: Mapping & QC (Steps 1-5)
- Calculate read statistics (N50, Q50, total bases)
- Map reads to reference genome with minimap2
- Filter alignments (MAPQ >= 20, remove secondary/supplementary)
- Generate alignment statistics and histograms

### Phase 2: Per-base Error (Steps 6-8)
- Calculate per-base error rates from filtered alignments (strand-specific)
- Split files by chromosome (checkpoint for parallelization)
- Correlate per-base error between samples via random subsampling

### Phase 3: Reactivity (Step 9)
- Calculate reactivity: treatment_error - control_error
- Operates per-chromosome for parallel processing

### Phase 4: Output Formats (Steps 10-14)
- Convert reactivity to bedGraph format with significance thresholds
- Calculate reactive nucleotide density in genomic windows
- Convert bedGraph to bigWig format
- Merge chromosome-split files back together

### Phase 5: Feature Annotation (Steps 15-18)
- Split feature BED files by chromosome
- Calculate signal around genomic features (TSS, etc.)
- Merge and average annotations

## Dependencies

Install via conda:
```bash
conda env create -f environment.yml
conda activate nanoprint
```

- samtools, minimap2, seqtk (alignment)
- bedtools, UCSC tools (bedGraphToBigWig, bigWigToBedGraph)
- gawk, bc, python3
- R with ggplot2, dplyr

## Directory Structure

```
CONFIG                   # User configuration file
Snakefile                # Auto-generated by CONFIG.sh
data/                    # Pipeline outputs
  aligned_reads/         # Raw mapped reads
  filtered_alignments/   # Quality-filtered BAMs
  perbase_error/         # Per-base error calculations
  perbase_error_by_chr/  # Split by chromosome
  reactivity/            # Treatment - control
  bg/                    # bedGraph files
  bw/                    # bigWig files
  bw_merged/             # Merged bigWig files
  windows/               # Density calculations
  annotations/           # Feature annotations
raw_data/                # Input FASTQ/BAM files
resources/
  genomes/               # Reference genomes (.fa)
  features/              # Feature BED files
tables/                  # Summary statistics
workflow/
  scripts/               # Shell scripts including CONFIG.sh
  rules/                 # Snakemake rule files (phases 1-5)
```

## Usage

1. Create CONFIG file with your samples and settings
2. Generate Snakefile: `./workflow/scripts/CONFIG.sh -i CONFIG -o Snakefile`
3. Run pipeline: `snakemake --cores <N>`

Or use individual scripts directly (see README for documentation).

---

## Development Notes

### Checkpoint Removal (Jan 2026)

The original workflow used Snakemake checkpoints (`split_perbase_by_chr`, `split_features_by_chr`) to dynamically determine chromosome wildcards at runtime. This caused issues:
- Could not generate DAG diagram before running
- `MissingInputException` errors when checkpoint outputs didn't match expected inputs
- Fragile dependency resolution

**Solution: Marker file + direct dependency pattern**

Instead of checkpoints, chromosome wildcards are now hardcoded at CONFIG.sh generation time:

1. **CONFIG.sh extracts chromosomes**:
   - `CHROMOSOMES` dict: from genome `.fai` files (`cut -f1`)
   - `FEATURE_CHROMOSOMES` dict: from BED column 1 (`cut -f1 | sort -u`), filtered to only chromosomes present in the genome

2. **Split rules output marker files**:
   ```python
   rule split_features_by_chr:
       output:
           done="resources/features/{feature}_by_chr/.done"
       shell:
           """
           python3 workflow/scripts/Split_by_chr.sh ...
           touch {output.done}
           """
   ```

3. **Downstream rules depend on `.done` marker, reference split files via `params`**:
   ```python
   rule annotate_features:
       input:
           feature_done="resources/features/{feature}_by_chr/.done",
           ...
       params:
           bed="resources/features/{feature}_by_chr/{feature}_{chr}.bed.gz"
       shell:
           "annotate_features.sh -b {params.bed} ..."
   ```

**Important**: Do NOT use "passthrough rules" that declare split files as `output` and just run `test -f`. Snakemake 8.x deletes output files before executing a rule, so the `test -f` will always fail. Instead, split files must be referenced via `params` (not tracked by Snakemake) with `.done` as the real dependency.

For phase 2 (perbase error), per-genome split rules in `genome_specific_rules.smk` declare all chromosome files as direct outputs, avoiding this issue entirely.

### Bash 3.2 Compatibility

macOS ships with bash 3.2 (due to GPL licensing). CONFIG.sh must avoid:
- **Associative arrays** (`declare -A`) - use parallel arrays instead
- **`mapfile`** - use `$(command)` with word splitting instead

Example workaround:
```bash
# Instead of associative arrays:
declare -a GENOME_CHR_NAMES=()
declare -a GENOME_CHR_VALUES=()
GENOME_CHR_NAMES+=("$genome")
GENOME_CHR_VALUES+=("$chrs")
```

### BED File Format

Feature BED files must be **BED6 format** with strand in **column 6** (not column 4):
```
chrom  start  end  name  score  strand
[0]    [1]    [2]  [3]   [4]    [5]
```

The `annotate_features.sh` script reads strand from `fields[5]`.

### File Path Conventions

Per-base error split files follow this structure:
```
data/perbase_error_by_chr/{genome}/{raw_sample}_{strand}/{raw_sample}_{strand}_{chr}.txt
```

BedGraph files (phase 4) are at:
```
data/bg/{genome}/{sample}_{strand}_{chr}.bg
```
(Note: NOT `data/bg_by_chr/` - this was a bug that was fixed)

### Chunk-Based Memory Management in annotate_features.sh (Feb 2026)

The `annotate_features.sh` script previously loaded entire per-base error files into memory, causing excessive memory usage for large chromosomes (~750MB+ for human chr1).

**Solution: Two-pass chunk-based processing with sliding window cache**

The script now requires a genome index file (`-g <genome.fai>`) and uses 100kb chunks:

**Pass 1: Split inputs into bin temp files (single read through each file)**
- Reads chromosome name from first line of treatment file (column 1)
- Looks up chromosome size from `.fai` file
- Splits each input file into 100kb bins based on position:
  - `bin_idx = position // 100000`
  - No data duplication between bins
- Creates temp files: `tmp/{treatment,control,reactivity}/bin_N.txt`, `tmp/features/bin_N.bed`

**Pass 2: Process features with sliding cache**
- Features are processed bin-by-bin in sorted order
- `ChunkCache` class keeps max 2 adjacent bins loaded at once
- For each feature:
  1. Calculate window span needed: `[ref_pos - window_span, ref_pos + window_span]`
  2. Determine which bins are needed (max 2 for boundary cases)
  3. Cache loads/unloads bins as needed
  4. Process windows using merged data from loaded bins

**Memory bound**: ~600KB max (2 bins × 3 files × 100kb) vs ~750MB+ previously

**Snakemake rule change**: `phase5_annotate_features.smk` now passes `fai="resources/genomes/{genome}.fa.fai"` to the script via `-g {input.fai}`

```python
# Sliding window example:
# Feature at position 98,000 needs data from 88,000 to 108,000
# - Bin 0 covers [0, 100000) → needed
# - Bin 1 covers [100000, 200000) → needed
# Cache loads both bins, processes feature, continues
```

### Single-Pass Averaging in average_feature_annotation.sh (Feb 2026)

The `average_feature_annotation.sh` script originally used a chunking approach that created individual files for each (distance, sample, strand) combination, then processed each file separately. This caused two problems:

**Bug**: Negative distance values (e.g., `-9990`) created filenames starting with `-`, causing awk redirection errors:
```
fatal: cannot redirect to '.../chunks/-9990_Treatment_rev.txt': Operation not permitted
```

**Performance**: With ~80,000 unique keys (±10,000 distances × 2 samples × 2 strands) and 100M+ input records:
- Created thousands of temp files
- Spawned thousands of awk subprocesses in a bash loop
- Massive I/O overhead

**Solution: Single-pass awk with associative arrays**

The script now uses a single awk invocation that:
1. Reads all input data in one pass
2. Accumulates sums using composite keys: `distance SUBSEP sample SUBSEP strand`
3. Outputs averages at the end
4. Pipes directly to sort and gzip

```bash
gunzip -c "$INPUT" | awk '
BEGIN { FS = "\t"; OFS = "\t" }
NR == 1 { next }
{
    key = $1 SUBSEP $5 SUBSEP $6
    sum_cov[key] += $2
    sum_err[key] += $3
    count[key]++
    if ($4 != "") { sum_react[key] += $4; react_count[key]++ }
}
END {
    print "Distance\tCoverage\tPerbase_error\tReactivity\tSample\tStrand"
    for (key in count) { ... output averages ... }
}
' | { read -r header; echo "$header"; sort -t$'\t' -k1,1n -k5,5 -k6,6; } | gzip -c > "$OUTPUT"
```

**Performance impact:**
- No temp files (eliminates I/O bottleneck)
- Single process (eliminates subprocess overhead)
- Memory: ~3-5MB for 80,000 keys (trivial)
- Runtime: seconds/minutes vs potentially hours

### Special Marker Values in Reactivity Files

Phase 3 (`Calculate_reactivity.sh`) outputs special marker values for positions missing data:
- `999999` = position missing in control file
- `-999999` = position missing in treatment file

**How downstream steps handle these markers:**

| Phase | Rule | Script | Handling |
|-------|------|--------|----------|
| 4 | `reactivity_to_bedgraph` | `react_to_bg.sh` | **Explicit filter** (line 97): `awk '$4 != 999999 && $4 != -999999'` |
| 4 | `reactivity_density` | `react_dens.sh` | Implicit - input is already filtered bedGraph |
| 4 | `bedgraph_to_bigwig` | `bg_to_bw.sh` | Implicit - input is already filtered bedGraph |
| 4 | `merge_density` | `Merge_density.sh` | No filtering needed - concatenates filtered data |
| 4 | `merge_bigwig` | `Merge_bigwig.sh` | No filtering needed - concatenates filtered data |
| 5 | `annotate_features` | `annotate_features.sh` | **Explicit skip** in Python: `if react == 999999 or react == -999999: continue` |
| 5 | `merge_annotations` | `Merge_annotations.sh` | No filtering needed - concatenates |
| 5 | `average_annotations` | `average_feature_annotation.sh` | No filtering needed - averages filtered data |

**Key points:**
1. Two explicit filter points: `react_to_bg.sh` (line 97) and `annotate_features.sh` (embedded Python)
2. All downstream scripts receive pre-filtered data
3. Windows with no valid reactivity data output empty string for reactivity column

### Snakemake Benchmarking (Feb 2026)

All 18 rules (excluding the lightweight `feature_chr_file` passthrough) have `benchmark:` directives to track resource usage.

**Benchmark output format (TSV):**
- `s` - Wall clock time (seconds)
- `h:m:s` - Human-readable time
- `max_rss` - Maximum resident set size (memory in MB)
- `max_vms` - Maximum virtual memory size (MB)
- `io_in` / `io_out` - I/O read/write (MB)
- `mean_load` - Mean CPU load
- `cpu_time` - Total CPU time (seconds)

**Benchmark directory structure:**
```
benchmarks/
├── phase1/
│   ├── genome_faidx/{genome}.tsv
│   ├── read_stats/{raw_sample}.tsv
│   ├── map_reads/{genome}/{raw_sample}.tsv
│   ├── filter_alignments/{genome}/{raw_sample}.tsv
│   ├── alignment_stats/{genome}/{raw_sample}.tsv
│   └── histograms/{alignment}/{genome}/{raw_sample}.tsv
├── phase2/
│   ├── perbase_error/{genome}/{raw_sample}_{strand}.tsv
│   ├── split_perbase_by_chr/{genome}/{raw_sample}_{strand}.tsv
│   └── correlation/{genome}/{raw_sample_a}_vs_{raw_sample_b}_{strand}.tsv
├── phase3/
│   └── calculate_reactivity/{genome}/{sample}_{strand}_{chr}.tsv
├── phase4/
│   ├── reactivity_to_bedgraph/{genome}/{sample}_{strand}_{chr}.tsv
│   ├── reactivity_density/{genome}/{sample}_{strand}_{chr}_{size}_{sig}.tsv
│   ├── merge_density/{genome}/{sample}_{strand}_{size}_{sig}.tsv
│   ├── bedgraph_to_bigwig/{genome}/{sample}_{strand}_{chr}_{sig}.tsv
│   └── merge_bigwig/{genome}/{sample}_{strand}_{sig}.tsv
└── phase5/
    ├── split_features_by_chr/{feature}.tsv
    ├── annotate_features/{genome}/{feature}/{sample}_{strand}_{chr}.tsv
    ├── merge_annotations/{genome}/{feature}/{sample}_{strand}.tsv
    └── average_annotations/{genome}/{feature}/{sample}_{strand}.tsv
```

**Resource-intensive rules to monitor:**
- `map_reads` - minimap2 alignment (hours for large datasets)
- `perbase_error` - full BAM processing per strand
- `annotate_features` - feature annotation with chunk-based memory management

### Genome-Specific Rules Generation (Feb 2026)

The `split_perbase_by_chr` rule needs to declare all chromosome output files at DAG build time, but different genomes have different chromosome sets. Snakemake's output blocks are evaluated at parse time before wildcards are resolved, so `CHROMOSOMES[wildcards.genome]` cannot be used in output declarations (unlike input functions like `get_density_chr_files`).

**Problem**: A generic rule with `{genome}` wildcard cannot dynamically determine its outputs based on the genome.

**Solution: Per-genome rule generation**

CONFIG.sh now generates `genome_specific_rules.smk` in the base directory (not tracked by git) containing one `split_perbase_by_chr_{genome}` rule per genome:

```python
# Auto-generated in genome_specific_rules.smk
rule split_perbase_by_chr_chicken_v23:
    input:
        error="data/perbase_error/chicken.v23/{raw_sample}_{strand}.txt.gz"
    output:
        [wrap_output("perbase_error_by_chr", f) for f in
         expand("data/perbase_error_by_chr/chicken.v23/{{raw_sample}}_{{strand}}/{{raw_sample}}_{{strand}}_{chr}.txt.gz",
                chr=CHROMOSOMES["chicken.v23"])]
    # ... shell section with genome hardcoded
```

**Key points:**
1. Each genome gets its own rule with the genome name hardcoded (not a wildcard)
2. Rule names use underscores instead of dots: `split_perbase_by_chr_chicken_v23` (dots replaced with `_`)
3. The chromosome list is specific to each genome from CHROMOSOMES dict
4. `genome_specific_rules.smk` is included by the Snakefile and regenerated each time CONFIG.sh runs
5. The generic `split_perbase_by_chr` rule was removed from `phase2_perbase_error.smk`

**Why not other approaches:**
- **Input functions for outputs**: Not supported by Snakemake - outputs must be determinable at parse time
- **ALL_CHROMOSOMES union**: Fails if genomes have different chromosome sets (missing files cause errors)
- **Marker file + passthrough rule**: Adds complexity with extra rules that just verify files exist

### Split_by_chr.sh Rewritten in Python (Feb 2026)

The chromosome-splitting script was rewritten from bash/awk to pure Python after two classes of bugs:

**Bug 1: macOS BSD awk silently drops output files** during parallel Snakemake execution. Affects both pipe-based (`print | "gzip > file"`) and direct file redirection (`print > file`). Non-deterministic — different chromosomes fail on different runs.

**Bug 2: Snakemake 8.x deletes output files before executing a rule.** The original `feature_chr_file` passthrough rule declared split files as `output` and ran `test -f` to verify them. Snakemake deleted the file before running the test, so it always failed. This was the actual cause of the persistent `test -f` failures — not filesystem sync issues.

**Split_by_chr.sh solution:**

The script (`workflow/scripts/Split_by_chr.sh`) now has a `#!/usr/bin/env python3` shebang and uses only Python — no awk, no bash `gzip`, no shell loops. It writes gzipped output directly via `gzip.open()`.

```python
# Core logic — single pass, writes .gz files directly
opener = gzip.open if compressed else open
files = {}
with opener(input_path, "rt") as fh:
    for line in fh:
        chr_name = line.split("\t", 1)[0]
        if chr_name not in files:
            files[chr_name] = gzip.open(f"{out_prefix}{chr_name}.{ext}.gz", "wt")
        files[chr_name].write(line)
for f in files.values():
    f.flush()
    f.close()
```

**Snakemake rule solution:**

Removed the `feature_chr_file` passthrough rule entirely. `annotate_features` now depends on `.done` as input and references the `.bed.gz` file via `params` (not tracked by Snakemake, so it won't be deleted):

```python
rule annotate_features:
    input:
        feature_done="resources/features/{feature}_by_chr/.done",  # dependency
        ...
    params:
        bed="resources/features/{feature}_by_chr/{feature}_{chr}.bed.gz",  # not managed by Snakemake
```

**Key points:**
- All file I/O is Python — no awk or bash subprocesses involved
- `gzip.open("wt")` writes compressed output directly (no intermediate uncompressed files)
- Handles unsorted input (non-contiguous chromosome blocks) correctly
- Output directory is cleaned before splitting to prevent stale file interference
- Same CLI interface (`-i`, `-d`, `-h`) — Snakemake rules call with `python3 workflow/scripts/Split_by_chr.sh`
- **Never use passthrough rules** (declare file as output + `test -f`) in Snakemake 8.x — use `.done` marker + `params` instead

### Temporary File Management (Feb 2026)

The workflow generates many intermediate files that can consume significant disk space. A configurable temporary file system was implemented using Snakemake's `temp()` wrapper.

**CONFIG syntax:**

Add `^t` prefix lines to mark directories as temporary:
```
^t  data/aligned_reads
^t  data/bg
^t  data/bw
^t  data/windows
^t  data/annotations
^t  data/annotations_merged
```

**How it works:**

1. **CONFIG.sh parses `^t` lines** and generates:
   - `TEMP_DIRS` list in Snakefile
   - `TEMP_OUTPUTS` dict mapping output keys to boolean
   - `wrap_output(key, path)` helper function

2. **Rule files use `wrap_output()`** to conditionally wrap outputs:
   ```python
   output:
       bam=wrap_output("aligned_reads_bam", "data/aligned_reads/{genome}/{raw_sample}.bam")
   ```

3. **Snakemake's `temp()`** automatically deletes files after all downstream rules complete

**Configurable behavior:**
- Adding a `^t` line marks that directory's files as temporary
- Removing the `^t` line keeps those files (no deletion)
- No `^t` lines = all intermediate files kept (original behavior)

**Output keys mapped to directories:**

| Key | Directory | Rule(s) |
|-----|-----------|---------|
| `aligned_reads_bam` | `data/aligned_reads` | `map_reads` |
| `perbase_error_by_chr` | `data/perbase_error_by_chr` | `split_perbase_by_chr` |
| `reactivity` | `data/reactivity` | `calculate_reactivity` |
| `bg` | `data/bg` | `reactivity_to_bedgraph` |
| `windows` | `data/windows` | `reactivity_density` |
| `bw` | `data/bw` | `bedgraph_to_bigwig` |
| `annotations` | `data/annotations` | `annotate_features` |
| `annotations_merged` | `data/annotations_merged` | `merge_annotations` |

**Protected files (never temporary):**
- `data/filtered_alignments/` - Quality-filtered BAMs (user-requested)
- `data/perbase_error/` - Final per-base error (non-split)
- `data/bw_merged/` - Final merged bigWig
- `data/windows_merged/` - Final merged density
- `data/annotations_averaged/` - Final averaged annotations
- `results/igv/` - IGV export files (user-facing outputs)
- `tables/` - All QC tables
- `logs/` and `benchmarks/` - Always preserved

### Phase 6: IGV Export (Feb 2026)

Phase 6 generates IGV-ready files from aligned BAMs. Controlled by two boolean flags in CONFIG:
- `^igv-bam`: Strand-split BAMs + indices
- `^igv-bigwig`: Coverage bigWig files for each strand-split BAM

**Variable-length CONFIG prefixes:**

The original CONFIG.sh parser extracted prefixes as exactly 2 characters (`${line:0:2}`). This was changed to `${line%%[[:space:]]*}` to support multi-character prefixes like `^igv-bam`. The values extraction was updated correspondingly from `${line:2}` to `${line#"$prefix"}`. This is fully backward-compatible with existing 2-character prefixes.

**Rules (in `workflow/rules/phase6_igv.smk`):**

| Rule | Description | Key command |
|------|-------------|-------------|
| `igv_split_bam` | Split BAM by strand | `samtools view -b -h {-F 0x10|-f 0x10}` |
| `igv_index_bam` | Create BAM index | `samtools index` |
| `igv_coverage_bigwig` | Coverage bigWig | `bedtools genomecov -ibam -bg` + `bedGraphToBigWig` |

**IGV source mapping:**

```python
IGV_SOURCE_DIRS = {
    "all_alignments": "data/aligned_reads",
    "filtered_alignments": "data/filtered_alignments",
}
```

The `{igv_source}` wildcard maps to different input directories. An input function `get_igv_source_bam()` resolves the source BAM path at runtime.

**Output structure:**
```
results/igv/
├── all_alignments/{genome}/{raw_sample}_{strand}.bam[.bai|.bw]
└── filtered_alignments/{genome}/{raw_sample}_{strand}.bam[.bai|.bw]
```

**Conditional targets in rule all:**

IGV targets use inline conditional lists:
```python
expand(...) if IGV_BAM else [],
expand(...) if IGV_BIGWIG else [],
```

These are generated by CONFIG.sh based on whether `^igv-bam` / `^igv-bigwig` flags are present.

**Empty BAM handling:**

`igv_coverage_bigwig` checks for reads with `samtools view | head -1 | grep -q .` before generating coverage. If the BAM is empty (no reads mapped to that genome for that strand), it creates an empty placeholder file instead of failing on `bedGraphToBigWig`.
