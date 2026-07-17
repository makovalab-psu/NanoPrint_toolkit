# NanoPrint_toolkit

# Overview

NanoPrint_toolkit is a Snakemake-based pipeline for analyzing chemical footprinting data from long-read (Oxford Nanopore) sequencing. The pipeline processes raw sequencing reads through alignment, per-base error quantification, reactivity calculation (treatment minus control), and downstream analyses including significance-filtered bigWig generation, reactive nucleotide density calculations, and feature annotation. The workflow is optimized for parallelization through chromosome-level splitting of data.

When raw pod5 files are provided as input, the pipeline additionally runs Dorado basecalling and Uncalled4 DTW signal alignment to produce a parallel **per-base signal deviation** track. This track reports the mean difference between the expected pore model current and the observed nanopore signal at each genomic position, enabling comparison of chemical footprinting signal via two independent approaches (base-call errors vs. raw ion current deviation).

The toolkit is composed of a series of scripts found in workflow/scripts. A user can use the pipeline as intended with Snakemake, or use individual scripts as documented below.

---

# Quick Start: GPU Preprocessing (`nanoprint preprocess`)

When sequencing data is raw pod5 files, the first steps — Dorado basecalling and Uncalled4 pore-model signal alignment — require a GPU and are computationally intensive. If your GPU cluster and CPU cluster are separate machines (e.g., a local GPU server and Penn State Roar), you can run these steps independently using the `nanoprint` command-line tool, then transfer the resulting Uncalled4 BAM to the CPU cluster and continue the Snakemake pipeline there.

## PATH setup (one time)

```bash
export PATH="/path/to/Nanoprint_toolkit/bin:$PATH"
# Add to ~/.bashrc or ~/.zshrc to persist
```

## Usage

```bash
nanoprint preprocess \
    -i /path/to/pod5/file_directory \
    -g /path/to/genome.fa \
    -o Output_bam_file.bam \
    -p <threads>
```

| Flag | Required | Description |
|------|----------|-------------|
| `-i` | yes | Pod5 file or directory (searched recursively for `.pod5` files) |
| `-g` | yes | Reference genome FASTA |
| `-o` | yes | Output Uncalled4 BAM file |
| `-p` | yes | Number of threads / parallel processes |
| `-m` | no | Dorado basecalling model (default: `sup`; also accepts `hac`, `fast`, or a full model name) |
| `-T` | no | Directory for intermediate files (default: next to output; deleted on exit) |

## What it does

| Step | Tool | Result |
|------|------|--------|
| 1 | Dorado | pod5 → basecalled BAM (with move tables via `--emit-moves`) |
| 2 | minimap2 | basecalled BAM → aligned BAM (move tags preserved via `-T`/`-y`) |
| 3 | samtools | aligned BAM → filtered BAM (MAPQ ≥ 20, no secondary/supplementary) |
| 4 | Uncalled4 | filtered BAM split by pod5 source → per-pod5 Uncalled4 BAMs → merged and sorted → **Uncalled4 BAM** |

**Why step 4 splits by pod5:** If uncalled4 runs against all pod5 files at once on a coordinate-sorted BAM, it must seek randomly across every pod5 file to retrieve each read's raw signal — resulting in severe I/O bottlenecking (observed: 4.6% CPU utilization over 9 days on a 36-thread job). Instead, `nanoprint preprocess`:

1. Checks for duplicate pod5 basenames and exits with an error if any are found (duplicate basenames would cause silent read loss, since dorado's `fn:Z:` tag stores only the basename)
2. Sorts the filtered BAM by the `fn:Z:` tag (`samtools sort -t fn`) so reads from the same pod5 are contiguous — this allows the split pass to keep only a single pipe open at a time, avoiding OS file-descriptor limits
3. Makes a single streaming pass through the fn-sorted BAM, routing reads into per-pod5 BAMs one at a time
4. Runs uncalled4 on each pod5 independently — each job reads one file sequentially, eliminating random I/O

**Temporary disk usage** during step 4 peaks at approximately **2× the filtered BAM size**: one copy for the fn-sorted intermediate BAM (`fn_sorted.bam` in `-T` temp dir) plus the accumulating split BAMs (deleted progressively as each pod5 batch completes). Make sure the `-T` temp directory has enough space before starting.

Intermediate files are cleaned up automatically. The output is two files:
- `Output_bam_file.bam` — coordinate-sorted Uncalled4 BAM with embedded DTW tags
- `Output_bam_file.bam.bai` — BAM index

## Continuing the pipeline on ROAR

Transfer both output files to ROAR, then point the `^r` line in CONFIG at the Uncalled4 BAM as if it were a standard pre-aligned BAM. The pipeline will use the Uncalled4 BAM for per-base error (Phase 2) and signal deviation (Phase 2b) without re-running Phase 0.

```
# CONFIG on ROAR — use the transferred Uncalled4 BAM directly
^r  MySample  /path/to/Output_bam_file.bam  /path/to/Control_uncalled4.bam
```

## Dependencies

| Tool | Install |
|------|---------|
| `dorado` | Download binary from https://github.com/nanoporetech/dorado/releases |
| `samtools`, `minimap2` | `conda install -c bioconda samtools minimap2` |
| `uncalled4` | `pip install setuptools==69.5.1 && pip install uncalled4 && pip install "pod5" "pyarrow>=14,<20"` |

See the [Dependencies](#dependencies) section for full installation details, including the critical `pyarrow<20` pin that prevents a deadlock in uncalled4.

---

# Inputs

## Raw Sequencing Reads

### Description
Oxford Nanopore sequencing reads for treatment and control samples. Treatment samples are typically treated with a chemical probe (e.g., permanganate) while control samples are untreated.

The pipeline auto-detects the input format at runtime — no CONFIG flag required:
- **FASTQ/BAM mode**: standard basecalled reads; produces the standard per-base error track only
- **Pod5 mode**: raw signal files; triggers Phase 0 (Dorado basecalling + Uncalled4 signal alignment) and adds a parallel per-base signal deviation track alongside the standard error track

### Format
Standard FASTQ (gzipped), unaligned BAM, or raw Oxford Nanopore pod5 files.

### Location

Paths are specified directly in the `^r` line of the CONFIG file — no fixed directory layout is required:

```
# In CONFIG (^r  Sample  Treatment  Control):
^r  MySample  /absolute/path/to/treatment.bam   /absolute/path/to/control.bam
^r  MySample  /sequencer/run01/pod5_pass        /sequencer/run02/pod5_pass
```

Supported input formats and how to specify them:

| Format | Path to provide | Mode |
|--------|----------------|------|
| Gzipped FASTQ | `/path/to/{sample}.fastq.gz` | Standard error track |
| Unaligned BAM | `/path/to/{sample}.bam` | Standard error track |
| Pod5 directory | `/path/to/pod5/run/` (searched recursively) | Pod5 mode — adds signal track |
| Single pod5 file | `/path/to/{sample}.pod5` | Pod5 mode — adds signal track |

The raw_sample name used throughout the pipeline is derived from the path basename minus extension (e.g., `/runs/Sample01.bam` → `Sample01`). If two paths from different `^r` lines would produce the same name, CONFIG.sh exits with an error listing the conflicting paths.

---

## Reference Genome

### Description
Reference genome sequence for read alignment.

### Format
Standard FASTA format with accompanying FAI index (created automatically by the pipeline).

### Location
`resources/genomes/{genome}.fa`

---

## Feature BED Files

### Description
Genomic features for annotation analysis (e.g., transcription start sites, G-quadruplexes). The pipeline calculates signal profiles around these features.

### Format
Standard BED6 format (tab-delimited, no header):
| Column | Name | Description |
|--------|------|-------------|
| 1 | chrom | Chromosome name |
| 2 | chromStart | Start position (0-based) |
| 3 | chromEnd | End position |
| 4 | name | Feature name/identifier |
| 5 | score | Score value |
| 6 | strand | Strand (+ or -) |

### Example
```
chr19_MATERNAL	5	26	GQ:HUNTER=-1.71	64	-
chr19_MATERNAL	36	57	GQ:HUNTER=-1.71	64	-
chr19_MATERNAL	60	81	GQ:HUNTER=-1.71	64	-
chr19_MATERNAL	84	105	GQ:HUNTER=-1.71	64	-
chr19_MATERNAL	109	130	GQ:HUNTER=-1.71	64	-
chr19_MATERNAL	133	154	GQ:HUNTER=-1.71	64	-
chr19_MATERNAL	175	194	GQ:HUNTER=-1.89	67	-
chr19_MATERNAL	197	218	GQ:HUNTER=-1.86	64	-
chr19_MATERNAL	221	242	GQ:HUNTER=-1.71	64	-
chr19_MATERNAL	245	266	GQ:HUNTER=-1.71	64	-
```

### Location
`resources/features/{feature}.bed`

---

# Outputs

## Per-base Error Files

### Description
Per-base error rates calculated from aligned reads, separated by strand. Reports the mean error probability at each genomic position along with quantiles that capture per-read variability.

### Format
Tab-delimited, gzipped, no header:
| Column | Name | Description |
|--------|------|-------------|
| 1 | chrom | Chromosome name |
| 2 | position | Genomic position (1-based) |
| 3 | nucleotide | Reference nucleotide (A, C, G, T) |
| 4 | coverage | Read coverage at position |
| 5 | error | Mean per-base error probability |
| 6 | q25 | 0.25 quantile of per-read error probabilities (lower 50% CI bound) |
| 7 | q75 | 0.75 quantile of per-read error probabilities (upper 50% CI bound) |
| 8 | q025 | 0.025 quantile of per-read error probabilities (lower 95% CI bound) |
| 9 | q975 | 0.975 quantile of per-read error probabilities (upper 95% CI bound) |

### Example
```
chr19_MATERNAL	1	C	45	0.022222	0.010000	0.031623	0.001585	0.063096
chr19_MATERNAL	2	C	47	0.021277	0.010000	0.031623	0.001585	0.063096
chr19_MATERNAL	3	T	48	0.020833	0.010000	0.031623	0.001585	0.063096
chr19_MATERNAL	4	A	52	0.019231	0.010000	0.025119	0.001585	0.063096
chr19_MATERNAL	5	A	55	0.018182	0.010000	0.025119	0.001585	0.050119
chr19_MATERNAL	6	C	58	0.017241	0.010000	0.025119	0.001585	0.050119
chr19_MATERNAL	7	C	60	0.016667	0.010000	0.025119	0.001585	0.050119
chr19_MATERNAL	8	C	62	0.016129	0.010000	0.025119	0.001585	0.050119
chr19_MATERNAL	9	T	65	0.015385	0.010000	0.025119	0.001585	0.050119
chr19_MATERNAL	10	A	67	0.014925	0.010000	0.025119	0.001585	0.050119
```

### Location
`data/perbase_error/{genome}/{sample}_{strand}.txt.gz`

---

## Per-base Signal Deviation Files (pod5 mode only)

### Description
Per-base pore model signal deviation computed from Uncalled4 DTW alignment. Available only when the raw input is pod5 files. Reports the mean `dtw.model_diff` (model current − observed current, in pA) across all reads at each genomic position, along with quantiles that capture per-read variability. Positive values indicate the observed ion current is lower than the pore model expectation.

### Format
Tab-delimited, gzipped, no header (10 columns):
| Column | Name | Description |
|--------|------|-------------|
| 1 | chrom | Chromosome name |
| 2 | position | Genomic position (1-based) |
| 3 | nucleotide | Reference nucleotide (A, C, G, T) |
| 4 | coverage | Number of reads contributing to this position |
| 5 | mean_deviation | Mean dtw.model_diff across reads (pA) |
| 6 | q25 | 0.25 quantile of per-read dtw.model_diff (lower 50% CI bound, pA) |
| 7 | q75 | 0.75 quantile of per-read dtw.model_diff (upper 50% CI bound, pA) |
| 8 | q025 | 0.025 quantile of per-read dtw.model_diff (lower 95% CI bound, pA) |
| 9 | q975 | 0.975 quantile of per-read dtw.model_diff (upper 95% CI bound, pA) |
| 10 | mean_sq | Mean squared deviation — mean(dtw.model_diff²) = sum(dtw.model_diff²)/N (pA²) |

### Location
`data/perbase_signal/{genome}/{sample}_{strand}.txt.gz`

---

## Signal Reactivity Files (pod5 mode only)

### Description
Signal reactivity calculated as the difference in mean squared pore model signal deviation between treatment and control at each position:

> signal_reactivity = mean(dtw.model_diff²)_treatment − mean(dtw.model_diff²)_control  (pA²)

Mean squared deviation (column 10 of perbase_signal files) captures the magnitude of signal perturbation regardless of sign, making it a more sensitive metric for detecting chemical modification-induced changes in ion current. Positive values indicate greater signal variance in the treatment relative to control.

### Format
Identical to reactivity files (4 columns: chrom, position, nucleotide, signal_reactivity).

### Location
`data/signal_reactivity/{genome}/{sample}_{strand}_{chr}.txt.gz`

---

## Signal BigWig Files (pod5 mode only)

### Description
Signal reactivity data in bigWig format. Filtered by significance threshold; merged across chromosomes.

### Format
UCSC bigWig binary format.

### Location
`data/signal_bw_merged/{genome}/significance_threshold_{sig}/{sample}_{strand}.bw`

---

## Reactivity Files

### Description
Reactivity values calculated as treatment error minus control error at each position. Positive values indicate increased error in treatment (chemical modification signal).

### Format
Tab-delimited, gzipped, no header:
| Column | Name | Description |
|--------|------|-------------|
| 1 | chrom | Chromosome name |
| 2 | position | Genomic position (1-based) |
| 3 | nucleotide | Reference nucleotide |
| 4 | reactivity | Treatment error - Control error (999999 = missing in control, -999999 = missing in treatment) |

### Example
```
chr19_MATERNAL	1	C	0.003241
chr19_MATERNAL	2	C	0.001852
chr19_MATERNAL	3	T	-0.000421
chr19_MATERNAL	4	A	0.002156
chr19_MATERNAL	5	A	0.004523
chr19_MATERNAL	6	C	0.001234
chr19_MATERNAL	7	C	-0.000156
chr19_MATERNAL	8	C	0.002789
chr19_MATERNAL	9	T	0.003456
chr19_MATERNAL	10	A	0.001567
```

### Location
`data/reactivity/{genome}/{sample}_{strand}_{chr}.txt.gz`

---

## BedGraph Files

### Description
Reactivity data in UCSC bedGraph format with significance threshold information in the header.

### Format
UCSC bedGraph (tab-delimited, with header comments):
| Column | Name | Description |
|--------|------|-------------|
| 1 | chrom | Chromosome name |
| 2 | chromStart | Start position (0-based) |
| 3 | chromEnd | End position |
| 4 | value | Reactivity value |

### Example
```
# Significance thresholds:
# p < 0.05 (black):    reactivity >= 0.002887
# p < 0.01 (#FF8C00):  reactivity >= 0.004327
# p < 0.001 (red):     reactivity >= 0.005892
# p < 1e-04 (#810000): reactivity >= 0.007456
chr19_MATERNAL	0	1	0.003241
chr19_MATERNAL	1	2	0.001852
chr19_MATERNAL	2	3	-0.000421
chr19_MATERNAL	3	4	0.002156
chr19_MATERNAL	4	5	0.004523
chr19_MATERNAL	5	6	0.001234
chr19_MATERNAL	6	7	-0.000156
chr19_MATERNAL	7	8	0.002789
chr19_MATERNAL	8	9	0.003456
chr19_MATERNAL	9	10	0.001567
```

### Location
`data/bg/{genome}/{sample}_{strand}_{chr}.bg`

---

## BigWig Files

### Description
Binary indexed format for efficient visualization of reactivity data in genome browsers. Filtered by significance threshold.

### Format
UCSC bigWig binary format (viewable in IGV, UCSC Genome Browser, etc.)

### Location
`data/bw_merged/{genome}/significance_threshold_{sig}/{sample}_{strand}.bw`

---

## IGV Strand-Split BAM Files

### Description
BAM files split by strand (forward and reverse) for IGV visualization. Generated from both raw aligned reads and quality-filtered alignments. Includes BAM indices (.bai) for random access. Enabled by adding `^igv-bam` to CONFIG.

### Format
Standard BAM format with BAM index (.bai).

### Location
```
results/igv/all_alignments/{genome}/{raw_sample}_{strand}.bam
results/igv/all_alignments/{genome}/{raw_sample}_{strand}.bam.bai
results/igv/filtered_alignments/{genome}/{raw_sample}_{strand}.bam
results/igv/filtered_alignments/{genome}/{raw_sample}_{strand}.bam.bai
```

---

## IGV Coverage BigWig Files

### Description
Coverage depth bigWig files generated from strand-split BAMs for IGV visualization. Enabled by adding `^igv-bigwig` to CONFIG.

### Format
UCSC bigWig binary format (viewable in IGV, UCSC Genome Browser, etc.)

### Location
```
results/igv/all_alignments/{genome}/{raw_sample}_{strand}.bw
results/igv/filtered_alignments/{genome}/{raw_sample}_{strand}.bw
```

---

## Mean Reactivity BedGraph/BigWig Files

### Description
Mean reactivity values averaged within genomic windows. Useful for visualizing smoothed reactivity signal in genome browsers.

### Format
Tab-delimited bedGraph (4 columns, no header):
| Column | Name | Description |
|--------|------|-------------|
| 1 | chrom | Chromosome name |
| 2 | chromStart | Window start (0-based) |
| 3 | chromEnd | Window end |
| 4 | mean_reactivity | Mean reactivity in window (0 if no data) |

### Example
```
chr19_MATERNAL	0	1000	0.001234
chr19_MATERNAL	1000	2000	0.002345
chr19_MATERNAL	2000	3000	-0.000123
chr19_MATERNAL	3000	4000	0.001567
chr19_MATERNAL	4000	5000	0.003456
```

### Location
- BedGraph: `data/bg_mean_merged/{genome}/window_size_{mean_size}/{sample}_{strand}.bg`
- BigWig: `data/bw_mean_merged/{genome}/window_size_{mean_size}/{sample}_{strand}.bw`

---

## Density BedGraph Files

### Description
Count of significantly reactive nucleotides within genomic windows.

### Format
Tab-delimited bedGraph, no header:
| Column | Name | Description |
|--------|------|-------------|
| 1 | chrom | Chromosome name |
| 2 | chromStart | Window start (0-based) |
| 3 | chromEnd | Window end |
| 4 | count | Number of reactive nucleotides in window |
| 5 | sum | Sum of reactivity values in window |

### Example
```
chr19_MATERNAL	0	1000000	523	2.456789
chr19_MATERNAL	1000000	2000000	612	3.123456
chr19_MATERNAL	2000000	3000000	489	2.234567
chr19_MATERNAL	3000000	4000000	534	2.567890
chr19_MATERNAL	4000000	5000000	601	2.890123
chr19_MATERNAL	5000000	6000000	578	2.678901
chr19_MATERNAL	6000000	7000000	545	2.345678
chr19_MATERNAL	7000000	8000000	567	2.456789
chr19_MATERNAL	8000000	9000000	589	2.567890
chr19_MATERNAL	9000000	10000000	612	2.789012
```

### Location
`data/windows_merged/{genome}/window_size_{size}/significance_threshold_{sig}/{sample}_{strand}.bg`

---

## Feature Annotation Files

### Description
Signal profiles around genomic features, showing coverage, per-base error, and reactivity as a function of distance from the feature.

### Format
Tab-delimited, gzipped, with header:
| Column | Name | Description |
|--------|------|-------------|
| 1 | Distance | Distance from feature reference point (bp) |
| 2 | Coverage | Average read coverage |
| 3 | Perbase_error | Average per-base error |
| 4 | Reactivity | Average reactivity (Treatment only) |
| 5 | Sample | Treatment or Control |
| 6 | Strand | for or rev |

### Example
```
Distance	Coverage	Perbase_error	Reactivity	Sample	Strand
-10000	45.23	0.018234	0.002341	Treatment	for
-9990	46.12	0.017892	0.002456	Treatment	for
-9980	44.89	0.018456	0.002234	Treatment	for
-9970	45.67	0.018123	0.002567	Treatment	for
-9960	46.34	0.017789	0.002678	Treatment	for
-10000	44.56	0.015893		Control	for
-9990	45.23	0.015456		Control	for
-9980	44.12	0.016012		Control	for
-9970	45.89	0.015556		Control	for
-9960	46.01	0.015111		Control	for
```

### Location
`data/annotations_averaged/{genome}/{feature}/{sample}_{strand}.txt.gz`

---

## Read Statistics Table

### Description
Aggregated read statistics from all samples in a single CSV table. Each row is one sample.

### Format
CSV with header:
| Column | Name | Description |
|--------|------|-------------|
| 1 | Sample | Sample name |
| 2 | Giga_bp | Total gigabases sequenced |
| 3 | Reads_million | Total reads in millions |
| 4 | N50 | Read length N50 (bp) |
| 5 | Q50 | Read quality Q50 (mean Phred at 50% of bases) |

### Example
```
Sample,Giga_bp,Reads_million,N50,Q50
Hsap_HG002_LCL_Mn04,12.34,0.82,18523,17
Hsap_HG002_LCL_CTRL,11.89,0.79,17891,16
```

### Location
`tables/read_stats_table.csv`

---

## Alignment Statistics Table

### Description
Aggregated alignment statistics from all samples in a single CSV table. Each sample produces two rows: one for raw (not filtered) alignments and one for quality-filtered alignments.

### Format
CSV with header:
| Column | Name | Description |
|--------|------|-------------|
| 1 | Sample | Sample name |
| 2 | Filter_status | `Not_filtered` (raw) or `Filtered` |
| 3 | Total_sequences | Number of reads |
| 4 | Total_length | Total bases aligned |
| 5 | Bases_mapped | Bases mapped to reference |
| 6 | Bases_mapped_cigar | Bases mapped (CIGAR-based) |
| 7 | Mismatches | Total mismatches |
| 8 | Error_rate | Per-base error rate |
| 9 | Average_length | Mean read length |
| 10 | Average_quality | Mean base quality |
| 11 | Primary_alignments | Primary alignment count |
| 12 | Secondary_alignments | Secondary alignment count |
| 13 | Supplementary_alignments | Supplementary alignment count |

### Location
`tables/alignment_stats_table.csv`

---

## Histogram PDFs

### Description
5-panel histogram PDFs showing alignment quality distributions for each sample and alignment type. Panels: Read Length (RL), Mapping Quality (MAPQ), Insertion size (INS), Deletion size (DEL), and Coverage (COV). Values are winsorized to 1st–99th percentile and re-binned into 50 uniform bins.

### Format
PDF (6×7 inches, 2-column × 3-row grid)

### Location
`plots/histograms/{alignment}/{genome}/{raw_sample}_histograms.pdf`

---

## Pairwise Correlation Table

### Description
Pairwise Spearman and Pearson correlations between all sample pairs. For each pair, both forward and reverse strand data are pooled before computing the correlation coefficients.

### Format
CSV with header:
| Column | Name | Description |
|--------|------|-------------|
| 1 | Sample_A | First sample name |
| 2 | Sample_B | Second sample name |
| 3 | Spearman | Spearman rank correlation (−1 to 1) |
| 4 | Pearson | Pearson correlation coefficient (−1 to 1) |

### Location
`tables/perbase_error_correlation/{genome}/Pairwise_correlation_table.csv`

---

## Pairwise Correlation Heatmap

### Description
Dual heatmap PDF displaying Spearman (left panel) and Pearson (right panel) pairwise correlation coefficients for all sample pairs. Values annotated in each tile. Color scale: turbo palette from −1 (blue) to 1 (red).

### Format
PDF (6×5 inches, two panels side-by-side)

### Location
`plots/perbase_error_correlation/{genome}/Pairwise_correlation_heatmap.pdf`

---

## Feature Annotation Plot

### Description
3-panel line plot PDF showing mean signal as a function of distance from genomic features. Panel 1: mean coverage. Panel 2: per-base error rate. Panel 3: reactivity. Lines are colored by Sample (Treatment=black, Control=grey) and styled by genome strand (forward=solid, reverse=dotted).

### Format
PDF (7×3 inches, 3 panels in a single row; legend on right of Reactivity panel)

### Location
`plots/annotations_averaged/{genome}/{feature}/{sample}.pdf`

---

# Dependencies

## Standard Installation

All standard dependencies can be installed via conda:

```bash
conda env create -f environment.yml
conda activate nanoprint
```

| Category | Tools |
|----------|-------|
| Workflow | snakemake |
| Alignment | samtools, minimap2, seqtk |
| Genomic intervals | bedtools, ucsc-bedgraphtobigwig, ucsc-bigwigtobedgraph |
| Text processing | gawk, bc |
| Languages | python (>=3.8), pandas, pysam, R |
| R packages | ggplot2, dplyr, gridExtra |

## Pod5 / Signal Analysis (Optional)

Required only when raw pod5 input files are used. These tools are **not** installable via conda and require separate setup.

### 1. Dorado Basecaller

Download a pre-built binary from the [Dorado releases page](https://github.com/nanoporetech/dorado/releases) and place it in your `PATH`:

```bash
# Example (Linux, v0.9):
wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.9.0-linux-x64.tar.gz
tar -xzf dorado-0.9.0-linux-x64.tar.gz
export PATH="$PWD/dorado-0.9.0-linux-x64/bin:$PATH"
```

GPU (CUDA) is used automatically if available. CPU basecalling is supported but slow.

### 2. Uncalled4 Signal Aligner

Install via pip **in the following exact order** (with the nanoprint conda environment active):

```bash
conda activate nanoprint

# REQUIRED first — newer setuptools breaks the uncalled4 build
pip install setuptools==69.5.1

pip install uncalled4

# CRITICAL: pin pyarrow
# lib-pod5 >=0.3.33 + pyarrow >=20 deadlocks inside uncalled4's C extension.
# lib-pod5==0.3.10 does not exist on PyPI; the latest available version is fine.
# Pinning pyarrow<20 is sufficient to prevent the deadlock.
pip install "pod5" "pyarrow>=14,<20"
```

Verify the installation:

```bash
uncalled4 --version
python -c "import pod5; print(pod5.__version__)"
```

# Instructions

## 1. Organize data

Make directories for inputs

```bash
mkdir -p resources/genomes resources/features
```

Place reference genomes in `resources/genomes/` and feature BED files in `resources/features/`. Raw sequencing reads can be located anywhere on the filesystem — paths are specified directly in the CONFIG `^r` line.

## 2. Edit CONFIG file

Your CONFIG file should look something like this:

```
#Genomes, place in resources/genomes/
^g	test_genome.fa

#Features to annotate, place in resources/features/
#^f	bed file	n_windows	window_size
^f	g4Discovery.bed	1000	10

#Mean bedgraph window sizes
^a	1000

#Density bedgraph Window size
^w	1000000
^w	10000

#Significance threshold level, 0 is all data (no threshold),  1 is p <= 0.05,  2 is p <= 0.01,  3 is p <= 0.001,  4 is p <= 0.0001.
#Multiple levels can be specified (one per line). Density files are not produced for level 0.
^s	2

#Relationships between files
# Treatment and control fields are absolute (or relative) paths to the input file or directory.
# For FASTQ/BAM: provide the path to the file.
# For pod5: provide the path to the pod5 directory (searched recursively for .pod5 files).
# The raw_sample name is derived from the path basename minus extension.
#^    Sample          Treatment                                        Control
^r	Hsap_HG002_LCL	/path/to/reads/Hsap_HG002_LCL_Mn04.bam	/path/to/reads/Hsap_HG002_LCL_CTRL.bam

# Pod5 example:
#^r	Hsap_HG002_LCL	/sequencer/output/run_Mn04/pod5_pass	/sequencer/output/run_CTRL/pod5_pass

# Temporary intermediate files (deleted after downstream rules complete)
# Remove these lines to keep all intermediate files
^t	data/aligned_reads
^t	data/bg
^t	data/bw
^t	data/windows
^t	data/annotations
^t	data/annotations_merged

#Other settings
^igv-bam
^igv-bigwig

# Dorado model for pod5 basecalling (optional; default: sup)
# Use a shorthand ('sup', 'hac', 'fast') or a full model name
^dorado-model sup
```
The wildcard variables are assigned designated as:

^g The genome you want to map to
^f Features to annotate with optional parameters: bed_file, n_windows (default: 1000), window_size (default: 10)
^a Window size for mean reactivity bedGraph/bigWig files
^w The window files you want in the windows bed files
^s The significance threshold for identifying reactive nucleotides (0 = all data, no threshold; 1–4 = p-value cutoffs; multiple lines allowed; level 0 does not produce density files)
^r The relationship between sequencing reads. Treatment and control fields are absolute (or relative) paths to the input file or directory. The raw_sample name is derived from the path basename minus extension.
^t Directories containing temporary files (auto-deleted after use)
^igv-bam Generate strand-split BAMs and indices for IGV visualization (flag, no value)
^igv-bigwig Generate coverage bigWig files for each strand-split BAM (flag, no value)
^dorado-model Dorado basecalling model for pod5 input (default: sup); applies to all pod5 samples

If you want to try out alternative variables, just add another row.

## 3. Generate Snakefile

```bash
./workflow/scripts/CONFIG.sh
```

## 4. Visualize the workflow to check the Snakefile

```bash
snakemake --dag | dot -Tpdf > dag.pdf
```

## 5. Run the pipeline

### Direct execution with Snakemake

```bash
snakemake --cores 4
```

When using pod5 input, add `--resources gpu=1` to prevent multiple Dorado jobs from running simultaneously and competing for GPU memory:

```bash
snakemake --cores 8 --resources gpu=1
```

This serializes `dorado_basecall` jobs (one at a time) while allowing all other rules to run in parallel. Each job uses all available GPUs; running two jobs concurrently causes CUDA out-of-memory errors.

### SLURM cluster execution (Penn State Roar)

Generate a ready-to-submit SLURM batch script using `SLURM_CONFIG.sh`. This runs a
`snakemake --dry-run` to count pending jobs, estimates wall time from benchmark data,
and produces an sbatch script with the appropriate resource requests.

```bash
# Generate the batch script
./workflow/scripts/SLURM_CONFIG.sh --env nanoprint

# Or with custom settings
./workflow/scripts/SLURM_CONFIG.sh \
    --env nanoprint \
    -i Snakefile \
    --cores 16 \
    --alloc open \
    -o my_submit.sh \
    --slog my_slurm_logs

# Submit to the cluster
sbatch 20260217_submit_nanoprint.sh
```

**Arguments:**

| Flag | Default | Description |
|------|---------|-------------|
| `--env` | (required) | Conda environment name or path |
| `-i` | `./Snakefile` | Input Snakefile |
| `-o` | `YYYYMMDD_submit_nanoprint.sh` | Output batch script |
| `--alloc` | `open` | Allocation (`open` = free queue; any other value uses `sla-prio`) |
| `--cores` | `8` | Number of cores (8 GB memory per core on standard partition) |
| `--slog` | `YYYYMMDD_slurm_logs` | Directory for SLURM log files |

The script automatically selects the SLURM partition based on `--alloc`: `open` uses the
free `open` partition, while a paid allocation ID routes to `sla-prio` with `--account`.

# Commands

---

## 0a. dorado_basecall (pod5 mode)

**Description:** Basecall raw Oxford Nanopore pod5 files using Dorado, emitting move tables required by Uncalled4.

**Script:** `workflow/scripts/Dorado_basecall.sh`

**Inputs:**
- Pod5 directory or single pod5 file (absolute path from `^r` CONFIG line)

**Outputs:**
- `data/basecalled/{sample}.bam` (unsorted, unaligned; sorted/aligned in subsequent map_reads step)

**Dependencies:**
- dorado (separate binary download; see Dependencies section)
- samtools

**Documentation:**

```
Usage: Dorado_basecall.sh -i <pod5_dir_or_file> -o <output.bam> -m <model> [-t <threads>]

Basecall Oxford Nanopore pod5 files using Dorado.
Emits move tables (--emit-moves) required by Uncalled4 signal alignment.

Required arguments:
    -i    Input: pod5 directory or single pod5 file
    -o    Output BAM file
    -m    Dorado model ('sup', 'hac', 'fast', or full model name)

Optional arguments:
    -t    Number of threads (default: 1; GPU usage controlled by dorado itself)
    -h    Show this help message

Notes:
    - GPU is used automatically if available; set CUDA_VISIBLE_DEVICES to control
    - Output BAM is unsorted; the pipeline sorts it after alignment in map_reads
    - --recursive is added automatically when the input is a directory, so pod5
      files nested in sequencer output subdirectories (e.g. pod5_pass/) are found

Example:
    Dorado_basecall.sh -i /absolute/path/to/pod5/Sample01/ -o data/basecalled/Sample01.bam -m sup -t 4
```

---

## 0b. uncalled4_align (pod5 mode)

**Description:** Align raw nanopore signal to the pore model using Uncalled4 DTW alignment. Produces a BAM with embedded DTW tags (compact; used by uncalled4_convert_tsv). This is the computationally expensive step — run once per sample.

**Script:** `workflow/scripts/Uncalled4_align.sh`

**Inputs:**
- `data/filtered_alignments/{genome}/{sample}.bam` (sequence-level alignments from map_reads)
- Pod5 directory or file (absolute path from `^r` CONFIG line; searched recursively)
- `resources/genomes/{genome}.fa`

**Outputs:**
- `data/uncalled4/{genome}/{sample}.bam` (Uncalled4 BAM with DTW tags)
- `data/uncalled4/{genome}/{sample}.bam.bai` (BAM index)

**Dependencies:**
- uncalled4 (pip install; see Dependencies section)
- samtools

**Resource requirements (WGS, from js4007 benchmarks):**
- Wall time: 6–16 minutes per sample per BAM file
- Peak RSS: 23–75 GB (request `mem_mb=80000` on SLURM)
- CPU load: 560–1078% (effectively 6–11 cores)

**Documentation:**

```
Usage: Uncalled4_align.sh -i <filtered.bam> -p <pod5_dir> -g <genome.fa> -o <out.bam> [-t <threads>]

Align raw nanopore signals to the pore model reference using Uncalled4 (BAM output).
Input BAM must have sequence-level alignments (from minimap2) and a move table
(--emit-moves from Dorado) so Uncalled4 can trace each read's signal.

Required arguments:
    -i    Input sequence-aligned BAM (filtered_alignments/{genome}/{sample}.bam)
    -p    Pod5 file directory or single pod5 file
    -g    Reference genome FASTA
    -o    Output Uncalled4 BAM with DTW signal alignment

Optional arguments:
    -t    Number of parallel processes for Uncalled4 (default: 8)
    -h    Show this help message

Notes:
    - uncalled4 returns non-zero when any reads fail DTW (even if most succeed);
      the script catches this with || true and verifies the output is non-empty.
    - Command syntax: uncalled4 align --bam-in <bam> --ref <fa> --reads <pod5> -p <N> -o <out.bam>

Example:
    Uncalled4_align.sh -i data/filtered_alignments/genome/Sample01.bam \
        -p /absolute/path/to/pod5/run01/ -g resources/genomes/genome.fa \
        -o data/uncalled4/genome/Sample01.bam -t 8
```

---

## 0c. uncalled4_convert_tsv (pod5 mode)

**Description:** Convert an Uncalled4 BAM to a strand-specific DTW TSV without re-running signal alignment. Pre-filters to one strand via `samtools view`, then calls `uncalled4 convert`. Run once per strand (for and rev) after `uncalled4_align`.

**Script:** `workflow/scripts/Uncalled4_convert_tsv.sh`

**Inputs:**
- `data/uncalled4/{genome}/{sample}.bam` (Uncalled4 BAM from `uncalled4_align`)
- `data/uncalled4/{genome}/{sample}.bam.bai`

**Outputs:**
- `data/uncalled4_tsv/{genome}/{sample}_{strand}.tsv`

**Dependencies:**
- uncalled4
- samtools

**Documentation:**

```
Usage: Uncalled4_convert_tsv.sh -i <uncalled4.bam> -o <out.tsv> -s <for|rev> [-t <threads>]

Convert an Uncalled4 signal-alignment BAM to a strand-specific DTW TSV.
Pre-filters the BAM to the requested strand with samtools view, then calls
"uncalled4 convert" to extract DTW columns without re-running alignment.

Required arguments:
    -i    Input Uncalled4 BAM (data/uncalled4/{genome}/{sample}.bam)
    -o    Output TSV (dtw metrics per read × reference position)
    -s    Strand: for (forward, -F 0x10) or rev (reverse, -f 0x10)

Optional arguments:
    -t    Number of parallel processes for uncalled4 convert (default: 4)
    -h    Show this help message

TSV columns extracted:
    dtw.current       Normalized mean read signal current (pA)
    dtw.current_sd    Signal current standard deviation
    dtw.start         Signal sample start index in raw trace
    dtw.length        Number of signal samples spanning this position
    dtw.model_diff    Model current - observed current (pA); positive = observed < expected
    dtw.base          Reference base (letter or integer-encoded: 0=A,1=C,2=G,3=T)

Example:
    Uncalled4_convert_tsv.sh -i data/uncalled4/genome/Sample01.bam \
        -o data/uncalled4_tsv/genome/Sample01_for.tsv -s for -t 4
```

---

## 1. read_stats

**Description:** Calculate sequencing statistics from raw reads.

**Script:** `workflow/scripts/Read_stats.sh`

**Inputs:**
- Raw reads file (path from `^r` CONFIG line: `.fastq.gz` or `.bam`)

**Outputs:**
- `tables/read_stats/{sample}.txt`

**Dependencies:**
- samtools (for BAM input)
- gawk
- bc

**Documentation:**

```
Usage: Read_stats.sh -i <input.fastq.gz|input.bam> -o <output.txt> [-T tmpdir]

Calculate read statistics from FASTQ or BAM files.

Required arguments:
    -i    Input file (*.fastq.gz or *.bam)
    -o    Output file (tab-delimited)

Optional arguments:
    -T    Temporary directory (default: same directory as output)
    -h    Show this help message

Input formats:
    - Gzipped FASTQ (*.fastq.gz)
    - Unmapped BAM (*.bam)

Output format (tab-delimited with header):
    1. File_prefix    - Input filename without extension
    2. Giga_bp        - Total gigabases
    3. Reads_million  - Total reads in millions
    4. N50            - Read length N50 (length at 50% of total bases)
    5. Q50            - Read quality Q50 (mean Phred at 50% of total bases)

Example:
    Read_stats.sh -i sample.fastq.gz -o sample_stats.txt
    Read_stats.sh -i sample.bam -o sample_stats.txt
```

---

## 2. map_reads

**Description:** Map raw reads to reference genome using minimap2.

**Script:** `workflow/scripts/Map_reads.sh`

**Inputs:**
- Raw reads file (path from `^r` CONFIG line: `.fastq.gz`, `.bam`, or `data/basecalled/{sample}.bam` for pod5)
- `resources/genomes/{genome}.fa`

**Outputs:**
- `data/aligned_reads/{genome}/{sample}.bam`

**Dependencies:**
- minimap2
- samtools

**Documentation:**

```
Usage: Map_reads.sh -i <input_file> -o <output_file> -g <genome.fasta> [-T <temp_dir>]

Map raw reads to reference genome using minimap2 and sort with samtools.

Required arguments:
    -i    Input file (*.fastq.gz or *.bam)
    -o    Output BAM file (sorted)
    -g    Reference genome FASTA file

Optional arguments:
    -T    Temporary directory

Example:
    Map_reads.sh -i sample.fastq.gz -g reference.fa -o sample.bam
```

---

## 3. filter_alignments

**Description:** Filter alignments by mapping quality and remove secondary/supplementary.

**Script:** `workflow/scripts/Filter_alignments.sh`

**Inputs:**
- `data/aligned_reads/{genome}/{sample}.bam`

**Outputs:**
- `data/filtered_alignments/{genome}/{sample}.bam`

**Dependencies:**
- samtools

**Documentation:**

```
Usage: Filter_alignments.sh -i <input_bam> -o <output_bam>

Filter BAM alignments by mapping quality (MAPQ >= 20) and remove
secondary (-F 0x100) and supplementary (-F 0x800) alignments.

Required arguments:
    -i    Input BAM file
    -o    Output filtered BAM file

Example:
    Filter_alignments.sh -i sample.bam -o sample_filtered.bam
```

---

## 4. alignment_stats

**Description:** Generate alignment statistics from raw and filtered BAMs.

**Script:** `workflow/scripts/Alignment_stats.sh`

**Inputs:**
- `data/aligned_reads/{genome}/{sample}.bam`
- `data/filtered_alignments/{genome}/{sample}.bam`

**Outputs:**
- `data/aligned_reads/{genome}/{sample}_flagstats.txt`
- `data/aligned_reads/{genome}/{sample}_stats.txt`
- `data/filtered_alignments/{genome}/{sample}_flagstats.txt`
- `data/filtered_alignments/{genome}/{sample}_stats.txt`
- `tables/alignment_stats/{genome}/{sample}.txt`

**Dependencies:**
- samtools

**Documentation:**

```
Usage: Alignment_stats.sh -a <raw.bam> -f <filtered.bam> -o <output.txt>

Generate alignment statistics from raw and filtered BAM files using samtools.

Required arguments:
    -a    Raw alignment BAM file
    -f    Filtered alignment BAM file
    -o    Output table file (tab-delimited)

Optional arguments:
    -h    Show this help message

Outputs:
    For each BAM file (raw and filtered):
        <bam>.bai           - BAM index file
        <bam>_stats.txt     - samtools stats output
        <bam>_flagstats.txt - samtools flagstat output

    Output table (tab-delimited with header):
        Column 1: Sample              - Sample name from BAM filename
        Column 2: Statistic           - Description of the statistic
        Column 3: Raw_alignment       - Value from raw alignment
        Column 4: Filtered_alignment  - Value from filtered alignment

Example:
    Alignment_stats.sh -a sample_raw.bam -f sample_filtered.bam -o sample_alignment_stats.txt
```

---

## 5. histograms

**Description:** Extract histogram data from samtools stats output.

**Script:** `workflow/scripts/Make_histograms.sh`

**Inputs:**
- `data/{alignment}/{genome}/{sample}_stats.txt`

**Outputs:**
- `data/{alignment}/{genome}/{sample}_histograms.txt`

**Dependencies:**
- gawk

**Documentation:**

```
Usage: Make_histograms.sh -i <input_stats.txt> -o <output_histograms.txt>

Extract histogram data from samtools stats output files.

Required arguments:
    -i    Input stats file (from samtools stats / Alignment_stats.sh)
    -o    Output histogram file (tab-delimited)

Optional arguments:
    -h    Show this help message

Output format (tab-delimited with header):
    Column 1: Var    - Histogram type (RL, MAPQ, INS, DEL, COV)
    Column 2: Value  - Bin value (length, quality, size, coverage)
    Column 3: Count  - Count for that bin

Histogram types extracted:
    RL   - Read length distribution
    MAPQ - Mapping quality distribution
    INS  - Insertion size distribution
    DEL  - Deletion size distribution
    COV  - Coverage distribution

Example:
    Make_histograms.sh -i sample_stats.txt -o sample_histograms.txt
```

---

## 6. perbase_error

**Description:** Calculate per-base error rates from filtered alignments.

**Script:** `workflow/scripts/perbase_error.sh`

**Inputs:**
- `data/filtered_alignments/{genome}/{sample}.bam`
- `resources/genomes/{genome}.fa`

**Outputs:**
- `data/perbase_error/{genome}/{sample}_{strand}.txt.gz`

**Dependencies:**
- samtools
- python3

**Documentation:**

```
Usage: perbase_error.sh <--rev|--for> -i <input.bam> -g <genome.fasta> -o <output.txt.gz> [-T tmpdir]

Calculate per-base error rates from filtered BAM alignments.

Positional argument (required, must be first):
    --rev    Process reverse strand reads only (FLAG 0x10 set)
    --for    Process forward strand reads only (FLAG 0x10 not set)

Required arguments:
    -i    Input BAM file (filtered/merged alignments)
    -g    Reference genome FASTA file
    -o    Output file (gzipped, *.txt.gz)

Optional arguments:
    -T    Temporary directory (default: same directory as output)
    -h    Show this help message

Output format (tab-delimited, gzipped):
    Column 1: Chromosome name
    Column 2: Position in chromosome (1-based)
    Column 3: Nucleotide identity
    Column 4: Coverage
    Column 5: Mean per-base error probability
    Column 6: Q25  — 0.25 quantile (lower 50% CI bound)
    Column 7: Q75  — 0.75 quantile (upper 50% CI bound)
    Column 8: Q025 — 0.025 quantile (lower 95% CI bound)
    Column 9: Q975 — 0.975 quantile (upper 95% CI bound)

Example:
    perbase_error.sh --for -i sample_filtered.bam -g reference.fasta -o sample_forward.txt.gz
    perbase_error.sh --rev -i sample_filtered.bam -g reference.fasta -o sample_reverse.txt.gz
```

---

## 7. split_perbase_by_chr

**Description:** Split per-base error file by chromosome for parallelization. (CHECKPOINT)

**Script:** `workflow/scripts/Split_by_chr.sh`

**Inputs:**
- `data/perbase_error/{genome}/{sample}_{strand}.txt.gz`

**Outputs:**
- `data/perbase_error_by_chr/{genome}/{sample}_{strand}/` (directory with per-chromosome files)

**Dependencies:**
- gawk

**Documentation:**

```
Usage: Split_by_chr.sh -i <input_file>

Split a tab-delimited file by chromosome (first column).
Output files are created in the same directory as input.

Required arguments:
    -i    Input file (tab-delimited with chromosome in column 1)
          Supported formats: .txt, .txt.gz, .bed, .bed.gz, .bg, .bg.gz

Optional arguments:
    -h    Show this help message

Output:
    Creates <prefix>_<chr>.<ext> for each chromosome in the same directory.
    Prefix is derived from input filename (without extension).

Example:
    Split_by_chr.sh -i data/perbase_error/sample_for.txt.gz
    # Creates: data/perbase_error/sample_for_chr1.txt, sample_for_chr2.txt, ...
```

---

## 8. correlation

**Description:** Correlate per-base error between two samples via random subsampling.

**Script:** `workflow/scripts/Correlation.sh`

**Inputs:**
- `data/perbase_error/{genome}/{sample_a}_{strand}.txt.gz`
- `data/perbase_error/{genome}/{sample_b}_{strand}.txt.gz`
- `resources/genomes/{genome}.fa.fai`

**Outputs:**
- `tables/correlation/{genome}/{sample_a}_vs_{sample_b}_{strand}.txt`

**Dependencies:**
- gawk
- coreutils (sort, shuf)

**Documentation:**

```
Usage: Correlation.sh -a <perbase_error1.txt.gz> -b <perbase_error2.txt.gz> -g <genome.fa.fai> -o <output.txt> [-s subsamples] [-T tmpdir]

Correlate per-base error rates between two samples by random subsampling.

Required arguments:
    -a    Per-base error file 1 (gzipped)
    -b    Per-base error file 2 (gzipped)
    -g    Genome sizes file (samtools faidx .fai format)
    -o    Output correlation file (tab-delimited)

Optional arguments:
    -s    Number of subsamples (default: 10000)
    -T    Temporary directory (default: same directory as output)
    -h    Show this help message

Input format (per-base error files, gzipped):
    Column 1: Chromosome
    Column 2: Position (1-based)
    Column 3: Nucleotide
    Column 4: Coverage
    Column 5: Per-base error

Output format (tab-delimited with header):
    Column 1: Chromosome
    Column 2: Nucleotide
    Column 3: Perbase_error_1
    Column 4: Perbase_error_2
    Column 5: Coverage_1
    Column 6: Coverage_2

Example:
    Correlation.sh -a sample1.txt.gz -b sample2.txt.gz -g genome.fa.fai -o correlation.txt -s 10000
```

---

## 2b. perbase_signal_deviation (pod5 mode)

**Description:** Compute per-base pore model signal deviation from an Uncalled4 DTW TSV. Groups per-read DTW measurements by reference position and computes multiple statistics from `dtw.model_diff` (model − observed current, pA). The 10-column output includes mean deviation, quantiles, and mean squared deviation (pA²). The mean squared deviation (column 10) is used as the signal metric for Phase 3b reactivity calculation.

**Script:** `workflow/scripts/perbase_signal_deviation.py`

**Inputs:**
- `data/uncalled4_tsv/{genome}/{sample}_{strand}.tsv` (strand-specific TSV from `uncalled4_convert_tsv`)
- `resources/genomes/{genome}.fa` (fallback for nucleotide lookup if dtw.base absent)

**Outputs:**
- `data/perbase_signal/{genome}/{sample}_{strand}.txt.gz`

**Dependencies:**
- python3, pandas, pysam (fallback only)

**Documentation:**

```
Usage: perbase_signal_deviation.py -i <dtw.tsv> -g <genome.fa> -o <output.txt.gz>

Compute per-base pore model signal deviation from an Uncalled4 DTW TSV.

Required arguments:
    -i    Input Uncalled4 DTW TSV (strand-filtered; produced by uncalled4_convert_tsv)
    -g    Reference genome FASTA (fallback for nucleotide lookup if dtw.base absent)
    -o    Output file (.txt.gz)

Optional arguments:
    -c    Minimum coverage to emit a position (default: 1)

Input TSV columns (subset used):
    dtw.model_diff    Model current - observed current (pA); the signal deviation metric
    dtw.base          Reference base (letter or integer; handled automatically)
    ref / seq_name    Chromosome name (column name varies by uncalled4 version)
    pos / seq_pos     0-based reference position (converted to 1-based in output)

Output format (tab-delimited, gzipped):
    Column 1:  Chromosome name
    Column 2:  Position (1-based)
    Column 3:  Nucleotide (from dtw.base; pysam FASTA as fallback)
    Column 4:  Coverage (reads at this position)
    Column 5:  Mean signal deviation (mean dtw.model_diff, pA)
    Column 6:  Q25  — 0.25 quantile of dtw.model_diff (lower 50% CI bound, pA)
    Column 7:  Q75  — 0.75 quantile of dtw.model_diff (upper 50% CI bound, pA)
    Column 8:  Q025 — 0.025 quantile of dtw.model_diff (lower 95% CI bound, pA)
    Column 9:  Q975 — 0.975 quantile of dtw.model_diff (upper 95% CI bound, pA)
    Column 10: Mean squared deviation — mean(dtw.model_diff^2) = sum(dtw.model_diff^2)/N (pA^2)

Notes:
    - dtw.model_diff = model - observed (positive = observed current lower than expected)
    - Positions where DTW failed (marked '*' in TSV) are excluded (read as NaN)
    - Handles uncalled4 version differences in column naming automatically

Example:
    perbase_signal_deviation.py -i data/uncalled4_tsv/genome/Sample01_for.tsv \
        -g resources/genomes/genome.fa -o data/perbase_signal/genome/Sample01_for.txt.gz
```

---

## 9. calculate_reactivity

**Description:** Calculate reactivity by comparing treatment to control samples (per chromosome).

**Script:** `workflow/scripts/Calculate_reactivity.sh`

**Inputs:**
- `data/perbase_error_by_chr/{genome}/{treatment_sample}_{strand}/{treatment_sample}_{strand}_{chr}.txt`
- `data/perbase_error_by_chr/{genome}/{control_sample}_{strand}/{control_sample}_{strand}_{chr}.txt`

**Outputs:**
- `data/reactivity/{genome}/{sample}_{strand}_{chr}.txt.gz`

**Dependencies:**
- python3
- gawk

**Documentation:**

```
Usage: Calculate_reactivity.sh -p <MnO4.txt.gz> -m <CTRL.txt.gz> -o <output.txt.gz> [-c threshold] [-T tmpdir]

Calculate reactivity from perbase error (treatment minus control).

Required arguments:
    -p    Treatment file (MnO4/plus, gzipped)
    -m    Control file (CTRL/minus, gzipped)
    -o    Output file (gzipped)

Optional arguments:
    -c    Minimum coverage threshold (default: 10)
    -T    Temporary directory (default: same directory as output)
    -h    Show this help message

Input format (5 columns, tab-separated):
    1. Chromosome name
    2. Position (1-based)
    3. Nucleotide identity
    4. Coverage
    5. Perbase error

Output format (4 columns, tab-separated):
    1. Chromosome name
    2. Position (1-based)
    3. Nucleotide identity
    4. Reactivity (treatment_error - control_error)
       Special values:
         999999  = position missing in control
        -999999  = position missing in treatment

Notes:
    - Input files must contain exactly one chromosome each
    - Chromosome names must match between treatment and control files
    - Chromosome size is determined from the maximum position in input files

Example:
    Calculate_reactivity.sh -p MnO4_chr1.txt.gz -m CTRL_chr1.txt.gz -o react_chr1.txt.gz
    Calculate_reactivity.sh -p MnO4_chr1.txt.gz -m CTRL_chr1.txt.gz -o react_chr1.txt.gz -c 20 -T /tmp
```

---

## 10. reactivity_to_bedgraph

**Description:** Convert reactivity data to bedGraph format.

**Script:** `workflow/scripts/react_to_bg.sh`

**Inputs:**
- `data/reactivity/{genome}/{sample}_{strand}_{chr}.txt.gz`

**Outputs:**
- `data/bg/{genome}/{sample}_{strand}_{chr}.bg`

**Dependencies:**
- gawk
- R (for quantile calculation)
- coreutils (shuf, sort)

**Documentation:**

```
Usage: react_to_bg.sh -i <reactivity.txt.gz> -o <output.bg> [-T tmpdir]

Convert reactivity data to bedGraph format with significance-based coloring.

Required arguments:
    -i    Input reactivity file (gzipped)
    -o    Output bedGraph file

Optional arguments:
    -T    Temporary directory (default: same directory as output)
    -h    Show this help message

Input format (4 columns, tab-separated):
    1. Chromosome name
    2. Position (1-based)
    3. Nucleotide identity
    4. Reactivity

Output format: UCSC bedGraph (4 columns: chr, start, end, reactivity)
    Header contains significance thresholds for reference:
    Thresholds: ns, p<0.05, p<0.01, p<0.001, p<1e-04
    Colors:     grey, black, #FF8C00, red, #810000

Example:
    react_to_bg.sh -i reactivity_chr1.txt.gz -o reactivity_chr1.bg
```

---

## 11. reactivity_density

**Description:** Calculate reactive nucleotide density in genomic windows (per chromosome).

**Script:** `workflow/scripts/react_dens.sh`

**Inputs:**
- `data/bg_by_chr/{genome}/{sample}_{strand}/{sample}_{strand}_{chr}.bg`
- `resources/genomes/{genome}.fa.fai`

**Outputs:**
- `data/windows/{genome}/window_size_{size}/significance_threshold_{sig}/{sample}_{strand}_{chr}.bg`

**Dependencies:**
- gawk
- bedtools

**Documentation:**

```
Usage: react_dens.sh -i <input.bg> -o <output.bg> -g <chrom.sizes.fai> [-w window] [-p sig_level] [-T tmpdir]

Calculate density of permanganate reactive nucleotides in genomic windows.

Required arguments:
    -i    Input bedGraph file (from react_to_bg.sh, with significance header)
    -o    Output bedGraph file (counts per window)
    -g    Chromosome sizes file (FAI format from samtools faidx)

Optional arguments:
    -w    Window size in bp (default: 1000)
    -p    Significance level filter (default: 4); must be 1–4 (level 0 not supported)
              1 = p <= 0.05
              2 = p <= 0.01
              3 = p <= 0.001
              4 = p <= 0.0001
    -T    Temporary directory (default: same directory as output)
    -h    Show this help message

Input format: bedGraph with significance thresholds in header
    Header example:
    # p < 0.05 (black):    reactivity >= 0.002886645
    # p < 0.01 (#FF8C00):  reactivity >= 0.004327449
    ...

Output format: bedGraph (5 columns)
    1. Chromosome name
    2. Start (0-based)
    3. End (1-based)
    4. Count of reactive nucleotides in window
    5. Sum of significant reactivity signal in window (after filtering)

Note: Significance level 0 (all data) is not supported for density calculation.
When ^s 0 is set in CONFIG, density files are not produced for that level.

Example:
    react_dens.sh -i react_chr1.bg -o density_chr1.bg -g genome.fa.fai -w 1000 -p 4
```

---

## 12. bedgraph_to_bigwig

**Description:** Convert bedGraph to bigWig format (per chromosome).

**Script:** `workflow/scripts/bg_to_bw.sh`

**Inputs:**
- `data/bg_by_chr/{genome}/{sample}_{strand}/{sample}_{strand}_{chr}.bg`
- `resources/genomes/{genome}.fa.fai`

**Outputs:**
- `data/bw/{genome}/significance_threshold_{sig}/{sample}_{strand}_{chr}.bw`

**Dependencies:**
- gawk
- bedGraphToBigWig (UCSC tools)

**Documentation:**

```
Usage: bg_to_bw.sh -i <input.bg> -o <output.bw> -g <chrom.sizes.fai> [-p sig_level] [-T tmpdir]

Convert bedGraph to bigWig format with optional significance filtering.

Required arguments:
    -i    Input bedGraph file (from react_to_bg.sh, with significance header)
    -o    Output bigWig file
    -g    Chromosome sizes file (FAI format from samtools faidx)

Optional arguments:
    -p    Significance level filter (default: 4)
              0 = return all data (no filtering)
              1 = p <= 0.05
              2 = p <= 0.01
              3 = p <= 0.001
              4 = p <= 0.0001
    -T    Temporary directory (default: same directory as output)
    -h    Show this help message

Input format: bedGraph with significance thresholds in header
    Header example:
    # p < 0.05 (black):    reactivity >= 0.002886645
    # p < 0.01 (#FF8C00):  reactivity >= 0.004327449
    ...

Output format: UCSC bigWig binary format

Example:
    bg_to_bw.sh -i react_chr1.bg -o react_chr1.bw -g genome.fa.fai -p 4
    bg_to_bw.sh -i react_chr1.bg -o react_chr1.bw -g genome.fa.fai -p 0  # all data
```

---

## 13. merge_bigwig

**Description:** Merge chromosome-split bigWig files in genome order.

**Script:** `workflow/scripts/Merge_bigwig.sh`

**Inputs:**
- `data/bw/{genome}/significance_threshold_{sig}/{sample}_{strand}_{chr}.bw` (multiple)
- `resources/genomes/{genome}.fa.fai`

**Outputs:**
- `data/bw_merged/{genome}/significance_threshold_{sig}/{sample}_{strand}.bw`

**Dependencies:**
- bigWigToBedGraph (UCSC tools)
- bedGraphToBigWig (UCSC tools)
- coreutils (sort)

**Documentation:**

```
Usage: Merge_bigwig.sh -g <genome.fa.fai> -o <output.bw> [-T tmpdir] <input1.bw> <input2.bw> ...

Merge chromosome-split bigWig files in genome order.

Required arguments:
    -g    Genome sizes file (samtools faidx .fai format)
    -o    Output bigWig file

Optional arguments:
    -T    Temporary directory (default: same directory as output)
    -h    Show this help message

Positional arguments:
    Remaining arguments are input bigWig files to merge

The input files are merged in the order of chromosomes as they appear
in the genome sizes file (first column of .fai).

Dependencies:
    - bigWigToBedGraph (UCSC tools)
    - bedGraphToBigWig (UCSC tools)

Example:
    Merge_bigwig.sh -g genome.fa.fai -o merged.bw sample_chr1.bw sample_chr2.bw sample_chrX.bw
```

---

## 14. merge_density

**Description:** Merge chromosome-split density files in genome order.

**Script:** `workflow/scripts/Merge_density.sh`

**Inputs:**
- `data/windows/{genome}/window_size_{size}/significance_threshold_{sig}/{sample}_{strand}_{chr}.bg` (multiple)
- `resources/genomes/{genome}.fa.fai`

**Outputs:**
- `data/windows_merged/{genome}/window_size_{size}/significance_threshold_{sig}/{sample}_{strand}.bg`

**Dependencies:**
- bash 4+ (for associative arrays)

**Documentation:**

```
Usage: Merge_density.sh -g <genome.fa.fai> -o <output.bg> <input1.bg> <input2.bg> ...

Merge chromosome-split density/bedgraph files in genome order.

Required arguments:
    -g    Genome sizes file (samtools faidx .fai format)
    -o    Output file

Positional arguments:
    Remaining arguments are input files to merge

Optional arguments:
    -h    Show this help message

The input files are merged in the order of chromosomes as they appear
in the genome sizes file (first column of .fai).

Example:
    Merge_density.sh -g genome.fa.fai -o merged.bg sample_chr1.bg sample_chr2.bg sample_chrX.bg
```

---

## 15. mean_reactivity_bedgraph

**Description:** Calculate mean reactivity in genomic windows (per chromosome).

**Script:** `workflow/scripts/react_mean_bg.sh`

**Inputs:**
- `data/reactivity/{genome}/{sample}_{strand}_{chr}.txt.gz`
- `resources/genomes/{genome}.fa.fai`

**Outputs:**
- `data/bg_mean/{genome}/window_size_{mean_size}/{sample}_{strand}_{chr}.bg`

**Dependencies:**
- bedtools
- gawk

**Documentation:**

```
Usage: react_mean_bg.sh -i <reactivity.txt.gz> -o <output.bg> -g <genome.fa.fai> [-w window] [-T tmpdir]

Calculate mean reactivity in genomic windows.

Required arguments:
    -i    Input reactivity file (gzipped, from Calculate_reactivity.sh)
    -o    Output bedGraph file (mean reactivity per window)
    -g    Chromosome sizes file (FAI format from samtools faidx)

Optional arguments:
    -w    Window size in bp (default: 1000)
    -T    Temporary directory (default: same directory as output)
    -h    Show this help message

Input format (4 columns, tab-separated):
    1. Chromosome name
    2. Position (1-based)
    3. Nucleotide identity
    4. Reactivity

Output format: bedGraph (4 columns)
    1. Chromosome name
    2. Start (0-based)
    3. End (1-based)
    4. Mean reactivity in window

Example:
    react_mean_bg.sh -i reactivity_chr1.txt.gz -o mean_chr1.bg -g genome.fa.fai -w 1000
```

---

## 16. merge_mean_bedgraph

**Description:** Merge chromosome-split mean reactivity bedGraph files in genome order.

**Script:** `workflow/scripts/Merge_density.sh` (reused)

**Inputs:**
- `data/bg_mean/{genome}/window_size_{mean_size}/{sample}_{strand}_{chr}.bg` (multiple)
- `resources/genomes/{genome}.fa.fai`

**Outputs:**
- `data/bg_mean_merged/{genome}/window_size_{mean_size}/{sample}_{strand}.bg`

**Dependencies:**
- coreutils (cat)

**Documentation:**

See [14. merge_density](#14-merge_density) — same script is used for chromosome-order concatenation.

---

## 17. mean_bedgraph_to_bigwig

**Description:** Convert merged mean reactivity bedGraph to bigWig format.

**Script:** `workflow/scripts/mean_bg_to_bw.sh`

**Inputs:**
- `data/bg_mean_merged/{genome}/window_size_{mean_size}/{sample}_{strand}.bg`
- `resources/genomes/{genome}.fa.fai`

**Outputs:**
- `data/bw_mean_merged/{genome}/window_size_{mean_size}/{sample}_{strand}.bw`

**Dependencies:**
- bedGraphToBigWig (UCSC tools)

**Documentation:**

```
Usage: mean_bg_to_bw.sh -i <input.bg> -o <output.bw> -g <genome.fa.fai> [-T tmpdir]

Convert a mean reactivity bedGraph file to bigWig format.

Required arguments:
    -i    Input bedGraph file (4-column: chr, start, end, mean_reactivity)
    -o    Output bigWig file
    -g    Chromosome sizes file (FAI format from samtools faidx)

Optional arguments:
    -T    Temporary directory (default: same directory as output)
    -h    Show this help message

Dependencies:
    - bedGraphToBigWig (UCSC tools)

Example:
    mean_bg_to_bw.sh -i mean_merged.bg -o mean_merged.bw -g genome.fa.fai
```

---

## 18. split_features_by_chr (Phase 5)

**Description:** Split feature BED file by chromosome for parallelization. (CHECKPOINT)

**Script:** `workflow/scripts/Split_by_chr.sh`

**Inputs:**
- `resources/features/{feature}.bed`

**Outputs:**
- `resources/features/{feature}_by_chr/` (directory with per-chromosome files)

**Dependencies:**
- gawk

**Documentation:**

```
Usage: Split_by_chr.sh -i <input_file>

Split a tab-delimited file by chromosome (first column).
Output files are created in the same directory as input.

Required arguments:
    -i    Input file (tab-delimited with chromosome in column 1)
          Supported formats: .txt, .txt.gz, .bed, .bed.gz, .bg, .bg.gz

Optional arguments:
    -h    Show this help message

Output:
    Creates <prefix>_<chr>.<ext> for each chromosome in the same directory.
    Prefix is derived from input filename (without extension).

Example:
    Split_by_chr.sh -i resources/features/TSS.bed
    # Creates: resources/features/TSS_chr1.bed, TSS_chr2.bed, ...
```

---

## 19. annotate_features

**Description:** Calculate signal around genomic features (per chromosome).

**Script:** `workflow/scripts/annotate_features.sh`

**Inputs:**
- `resources/features/{feature}_by_chr/{feature}_{chr}.bed`
- `data/perbase_error_by_chr/{genome}/{treatment_sample}_{strand}/{treatment_sample}_{strand}_{chr}.txt`
- `data/perbase_error_by_chr/{genome}/{control_sample}_{strand}/{control_sample}_{strand}_{chr}.txt`
- `data/reactivity/{genome}/{sample}_{strand}_{chr}.txt.gz`

**Outputs:**
- `data/annotations/{genome}/{feature}/{sample}_{strand}_{chr}.txt.gz`

**Dependencies:**
- python3
- gawk

**Documentation:**

```
Usage: annotate_features.sh -b <annotations.bed> -t <treatment.txt> -c <control.txt> -r <reactivity.txt> -o <output.txt.gz> [OPTIONS]

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
    annotate_features.sh -b features_chr1.bed -t treat_for_chr1.txt -c ctrl_for_chr1.txt -r react_for_chr1.txt -o output_chr1.txt.gz
```

---

## 20. merge_annotations

**Description:** Merge chromosome-split annotation files.

**Script:** `workflow/scripts/Merge_annotations.sh`

**Inputs:**
- `data/annotations/{genome}/{feature}/{sample}_{strand}_{chr}.txt.gz` (multiple)

**Outputs:**
- `data/annotations_merged/{genome}/{feature}/{sample}_{strand}.txt.gz`

**Dependencies:**
- gzip

**Documentation:**

```
Usage: Merge_annotations.sh -o <output.txt.gz> <input1.txt.gz> <input2.txt.gz> ...

Merge chromosome-split annotation files into a single file.

Required arguments:
    -o    Output file (gzipped)

Positional arguments:
    Remaining arguments are input annotation files (gzipped)

Optional arguments:
    -h    Show this help message

Example:
    Merge_annotations.sh -o merged.txt.gz chr1.txt.gz chr2.txt.gz chr3.txt.gz
```

---

## 21. average_annotations

**Description:** Average annotations by distance, sample, and strand.

**Script:** `workflow/scripts/average_feature_annotation.sh`

**Inputs:**
- `data/annotations_merged/{genome}/{feature}/{sample}_{strand}.txt.gz`

**Outputs:**
- `data/annotations_averaged/{genome}/{feature}/{sample}_{strand}.txt.gz`

**Dependencies:**
- gawk
- gzip

**Documentation:**

```
Usage: average_feature_annotation.sh -i <input.txt.gz> -o <output.txt.gz> [-T tmpdir]

Average feature annotations by distance, sample (Treatment/Control), and strand.
Uses chunking for memory-efficient processing of large files.

Required arguments:
    -i    Input merged annotation file (gzipped)
    -o    Output averaged file (gzipped)

Optional arguments:
    -T    Temporary directory (default: output directory)
    -h    Show this help message

Input format (tab-delimited, gzipped):
    Distance, Coverage, Perbase_error, Reactivity, Sample, Strand

Output format (tab-delimited, gzipped):
    Distance, Coverage, Perbase_error, Reactivity, Sample, Strand
    (averaged across all features for each Distance/Sample/Strand combination)

Example:
    average_feature_annotation.sh -i merged_annotations.txt.gz -o averaged_annotations.txt.gz
```

---

## 22. igv_split_bam

**Description:** Split BAM into forward and reverse strand reads for IGV visualization.

**Inputs:**
- `data/aligned_reads/{genome}/{raw_sample}.bam` or `data/filtered_alignments/{genome}/{raw_sample}.bam`

**Outputs:**
- `results/igv/{all_alignments,filtered_alignments}/{genome}/{raw_sample}_{strand}.bam`

**Dependencies:**
- samtools

**Usage:**
Enabled by adding `^igv-bam` to CONFIG. Uses `samtools view -b -h -F 0x10` (forward) and `samtools view -b -h -f 0x10` (reverse) to split reads by mapping strand.

---

## 23. igv_index_bam

**Description:** Create BAM index for IGV visualization.

**Inputs:**
- `results/igv/{all_alignments,filtered_alignments}/{genome}/{raw_sample}_{strand}.bam`

**Outputs:**
- `results/igv/{all_alignments,filtered_alignments}/{genome}/{raw_sample}_{strand}.bam.bai`

**Dependencies:**
- samtools

**Usage:**
Enabled by adding `^igv-bam` to CONFIG. Uses `samtools index`.

---

## 24. igv_coverage_bigwig

**Description:** Generate coverage bigWig from strand-split BAM for IGV visualization.

**Inputs:**
- `results/igv/{all_alignments,filtered_alignments}/{genome}/{raw_sample}_{strand}.bam`
- `results/igv/{all_alignments,filtered_alignments}/{genome}/{raw_sample}_{strand}.bam.bai`
- `resources/genomes/{genome}.fa.fai`

**Outputs:**
- `results/igv/{all_alignments,filtered_alignments}/{genome}/{raw_sample}_{strand}.bw`

**Dependencies:**
- bedtools
- bedGraphToBigWig (UCSC tools)

**Usage:**
Enabled by adding `^igv-bigwig` to CONFIG. Generates coverage bedGraph via `bedtools genomecov -ibam -bg`, then converts to bigWig with `bedGraphToBigWig`. Creates an empty file if the BAM has no reads.

---

## 25. read_stats_table

**Description:** Combine per-sample read statistics into a single CSV table.

**Script:** `workflow/scripts/read_stats_table.sh`

**Inputs:**
- `tables/read_stats/{raw_sample}.txt` (all samples)

**Outputs:**
- `tables/read_stats_table.csv`

**Dependencies:**
- gawk

**Documentation:**

```
Usage: read_stats_table.sh -o <output.csv> <input1.txt> [input2.txt ...]

Combine per-sample read statistics into a single CSV table.

Required arguments:
    -o    Output CSV file

Positional arguments:
    One or more read_stats txt files (tab-delimited, output of Read_stats.sh)

Output format (CSV with header):
    Sample, Giga_bp, Reads_million, N50, Q50

Example:
    read_stats_table.sh -o tables/read_stats_table.csv \
        tables/read_stats/sample1.txt tables/read_stats/sample2.txt
```

---

## 26. alignment_stats_table

**Description:** Combine per-sample alignment statistics into a single CSV table (two rows per sample: raw and filtered).

**Script:** `workflow/scripts/alignment_stats_table.sh`

**Inputs:**
- `tables/alignment_stats/{genome}/{raw_sample}.txt` (all samples)

**Outputs:**
- `tables/alignment_stats_table.csv`

**Dependencies:**
- gawk

**Documentation:**

```
Usage: alignment_stats_table.sh -o <output.csv> <input1.txt> [input2.txt ...]

Combine per-sample alignment statistics into a single CSV table.
Each input file produces two rows: one for raw alignments, one for filtered.

Required arguments:
    -o    Output CSV file

Positional arguments:
    One or more alignment_stats txt files (tab-delimited, output of Alignment_stats.sh)

Input format (tab-delimited):
    Sample  Statistic  Raw_alignment  Filtered_alignment

Output format (CSV with header):
    Sample, Filter_status, Total_sequences, Total_length, Bases_mapped,
    Bases_mapped_cigar, Mismatches, Error_rate, Average_length, Average_quality,
    Primary_alignments, Secondary_alignments, Supplementary_alignments

Filter_status values:
    Not_filtered  - Raw alignment statistics
    Filtered      - Quality-filtered alignment statistics

Example:
    alignment_stats_table.sh -o tables/alignment_stats_table.csv \
        tables/alignment_stats/genome1/sample1.txt \
        tables/alignment_stats/genome1/sample2.txt
```

---

## 27. summary_histogram_plot

**Description:** Generate 5-panel histogram PDF from samtools histogram data (Read Length, MAPQ, Insertion size, Deletion size, Coverage).

**Script:** `workflow/scripts/Plot_summary_histograms.R`

**Inputs:**
- `data/{alignment}/{genome}/{raw_sample}_histograms.txt`

**Outputs:**
- `plots/histograms/{alignment}/{genome}/{raw_sample}_histograms.pdf`

**Dependencies:**
- R with ggplot2, dplyr, gridExtra

**Documentation:**

```
Usage: Rscript Plot_summary_histograms.R <input_histograms.txt> <output.pdf>

Generate a 5-panel histogram PDF from alignment histogram data.

Arguments:
    1    Input histogram file (tab-delimited, from Make_histograms.sh)
    2    Output PDF file

Input format (tab-delimited with header):
    Var    - Histogram type (RL, MAPQ, INS, DEL, COV)
    Value  - Bin value
    Count  - Count for that bin

Output:
    PDF (6×7 inches) with 5 panels arranged in a 2-column × 3-row grid:
        RL   - Read Length distribution
        MAPQ - Mapping Quality distribution
        INS  - Insertion size distribution
        DEL  - Deletion size distribution
        COV  - Coverage distribution
    Each panel applies 99% winsorization and 50-bin re-binning.

Example:
    Rscript Plot_summary_histograms.R sample_histograms.txt sample_histograms.pdf
```

---

## 28. summarize_correlation

**Description:** Compute pairwise Spearman/Pearson correlations from all per-sample-pair correlation files and generate a dual heatmap PDF.

**Script:** `workflow/scripts/Summarize_correlation.R`

**Inputs:**
- `tables/correlation/{genome}/{raw_sample_a}_vs_{raw_sample_b}_{strand}.txt` (all pairs and strands)

**Outputs:**
- `tables/perbase_error_correlation/{genome}/Pairwise_correlation_table.csv`
- `plots/perbase_error_correlation/{genome}/Pairwise_correlation_heatmap.pdf`

**Dependencies:**
- R with ggplot2, dplyr, gridExtra

**Documentation:**

```
Usage: Rscript Summarize_correlation.R <out_table.csv> <out_plot.pdf> <input1.txt> [input2.txt ...]

Compute pairwise Spearman and Pearson correlations and generate heatmap plots.

Arguments:
    1    Output CSV table
    2    Output PDF heatmap
    3+   Input correlation files (output of Correlation.sh, any number)

Input format (tab-delimited with header):
    Chromosome, Nucleotide, Perbase_error_1, Perbase_error_2, Coverage_1, Coverage_2

Output table (CSV):
    Sample_A, Sample_B, Spearman, Pearson

Output plot (PDF, 6×5 inches):
    Two-panel heatmap: Spearman correlation (left) and Pearson correlation (right).
    Values annotated in each tile. Color scale: turbo palette from -1 to 1.

Notes:
    - Both forward and reverse strand files for each pair are pooled before computing
      correlation coefficients.
    - Sample names are parsed from filenames: {sample_a}_vs_{sample_b}_{strand}.txt

Example:
    Rscript Summarize_correlation.R \
        tables/perbase_error_correlation/genome/Pairwise_correlation_table.csv \
        plots/perbase_error_correlation/genome/Pairwise_correlation_heatmap.pdf \
        tables/correlation/genome/sample1_vs_sample2_for.txt \
        tables/correlation/genome/sample1_vs_sample2_rev.txt
```

---

## 29. plot_annotation

**Description:** Plot mean coverage, per-base error, and reactivity as a function of distance from genomic features.

**Script:** `workflow/scripts/Plot_annotation.R`

**Inputs:**
- `data/annotations_averaged/{genome}/{feature}/{sample}_for.txt.gz`
- `data/annotations_averaged/{genome}/{feature}/{sample}_rev.txt.gz`

**Outputs:**
- `plots/annotations_averaged/{genome}/{feature}/{sample}.pdf`

**Dependencies:**
- R with ggplot2, dplyr, gridExtra

**Documentation:**

```
Usage: Rscript Plot_annotation.R <forward.txt.gz> <reverse.txt.gz> <output.pdf>

Generate a 3-panel line plot PDF from averaged annotation files.

Arguments:
    1    Forward strand averaged annotation file (gzipped)
    2    Reverse strand averaged annotation file (gzipped)
    3    Output PDF file

Input format (tab-delimited, gzipped, with header):
    Distance, Coverage, Perbase_error, Reactivity, Sample, Strand

Output (PDF, 7×3 inches, 3 panels in a single row; legend on right of Reactivity panel):
    Panel 1: Mean Coverage vs distance
    Panel 2: Per-Base Error vs distance
    Panel 3: Reactivity vs distance

Aesthetics:
    Color:    Treatment = black, Control = grey
    Linetype: Forward strand = solid, Reverse strand = dotted

Notes:
    - The BED feature strand column (+ / -) is averaged away; only the genome
      strand (forward/reverse, from filename) affects line style.
    - Reactivity is only plotted for Treatment rows.

Example:
    Rscript Plot_annotation.R \
        data/annotations_averaged/genome/feature/sample_for.txt.gz \
        data/annotations_averaged/genome/feature/sample_rev.txt.gz \
        plots/annotations_averaged/genome/feature/sample.pdf
```

---

# Technical Notes

## Special Marker Values in Reactivity Files

Phase 3 (`Calculate_reactivity.sh`) outputs special marker values for positions where data is missing:

| Value | Meaning |
|-------|---------|
| `999999` | Position missing in control file |
| `-999999` | Position missing in treatment file |

These markers are filtered out at two points in the pipeline:

1. **Phase 4 - `react_to_bg.sh`**: Removes rows with marker values before converting to bedGraph format
2. **Phase 5 - `annotate_features.sh`**: Skips marker values when calculating window averages

All downstream steps (density calculation, bigWig conversion, merging, averaging) receive pre-filtered data and require no special handling. Windows that contain no valid reactivity data output an empty string for the reactivity column rather than a default value.

## Temporary File Management

The pipeline can automatically delete intermediate files after they are no longer needed, reducing disk usage. This is controlled via `^t` prefix lines in the CONFIG file.

### How It Works

Files in directories marked with `^t` are wrapped with Snakemake's `temp()` function, which automatically deletes them after all downstream rules that depend on them have completed.

### Temporary vs Kept Files

**Files marked as temporary (deleted after use):**

| Directory | Contents | Deleted After |
|-----------|----------|---------------|
| `data/aligned_reads/` | Raw BAM alignments | `filter_alignments`, `alignment_stats` complete |
| `data/bg/` | Per-chromosome bedGraph | `reactivity_density`, `bedgraph_to_bigwig` complete |
| `data/bw/` | Per-chromosome bigWig | `merge_bigwig` complete |
| `data/windows/` | Per-chromosome density | `merge_density` complete |
| `data/bg_mean/` | Per-chromosome mean reactivity bedGraph | `merge_mean_bedgraph` complete |
| `data/annotations/` | Per-chromosome annotations | `merge_annotations` complete |
| `data/annotations_merged/` | Pre-averaged annotations | `average_annotations` complete |

**Files always kept (final outputs):**

| Directory | Contents |
|-----------|----------|
| `data/filtered_alignments/` | Quality-filtered BAM files |
| `data/perbase_error/` | Per-base error rates (whole genome) |
| `data/perbase_error_by_chr/` | Per-base error (chromosome splits) |
| `data/reactivity/` | Reactivity values (chromosome splits) |
| `data/bw_merged/` | Final merged bigWig files |
| `data/bg_mean_merged/` | Final merged mean reactivity bedGraph |
| `data/bw_mean_merged/` | Final merged mean reactivity bigWig |
| `data/windows_merged/` | Final merged density files |
| `data/annotations_averaged/` | Final averaged annotations |
| `tables/` | All QC statistics tables |
| `logs/` | Execution logs |
| `benchmarks/` | Performance metrics |

### Customizing Temporary Files

To keep all intermediate files (original behavior), remove all `^t` lines from CONFIG.

To mark additional directories as temporary, add `^t` lines:
```
^t  data/perbase_error_by_chr
^t  data/reactivity
```

After modifying CONFIG, regenerate the Snakefile:
```bash
./workflow/scripts/CONFIG.sh -i CONFIG -o Snakefile
```

## Genome-Specific Rules

The pipeline supports multiple genomes with different chromosome sets. Because Snakemake requires all output files to be declared at parse time (before wildcards are resolved), rules that split files by chromosome cannot use a generic `{genome}` wildcard to determine the chromosome list dynamically.

### How It Works

When you run `CONFIG.sh`, it generates two files:
1. **Snakefile** - The main workflow file with configuration and includes
2. **genome_specific_rules.smk** - Auto-generated rules for chromosome splitting, one per genome

For example, with `chicken.v23` genome, CONFIG.sh generates:
```python
rule split_perbase_by_chr_chicken_v23:
    input:
        error="data/perbase_error/chicken.v23/{raw_sample}_{strand}.txt.gz"
    output:
        # Expands to all 41 chicken chromosomes
        [wrap_output("perbase_error_by_chr", f) for f in
         expand("data/perbase_error_by_chr/chicken.v23/..._{chr}.txt.gz",
                chr=["chr1", "chr2", ..., "chrZ", "mito"])]
```

### Important Notes

- **Always regenerate after CONFIG changes**: Run `./workflow/scripts/CONFIG.sh` whenever you modify the CONFIG file
- **genome_specific_rules.smk is not tracked by git**: This file is generated fresh for each implementation and should not be committed
- **Multi-genome support**: Each genome in CONFIG gets its own splitting rule with the correct chromosome list

## Pod5 Auto-Detection and Signal Analysis Track

### How Pod5 Mode Is Triggered

The pipeline detects pod5 input automatically at Snakemake run time via the `has_pod5()` function (defined in `workflow/rules/phase0_pod5_processing.smk`). The path for each raw sample comes from the `RAW_PATHS` dict generated by CONFIG.sh from the `^r` line:

1. If the path ends in `.pod5` and is a file → pod5 mode
2. If the path is a directory and any `.pod5` file exists anywhere under it (`os.walk`, recursive) → pod5 mode
3. Otherwise → standard FASTQ/BAM mode

No CONFIG flag is needed. The `find_raw_reads()` function in phase 1 returns `data/basecalled/{sample}.bam` when pod5 is detected, causing Snakemake to add `dorado_basecall` as an upstream dependency automatically.

### Signal Analysis Data Flow

The signal deviation track runs in parallel with the standard error track:

```
pod5 → dorado_basecall → (basecalled BAM)
                              ↓
                     map_reads + filter_alignments
                              ↓
              ┌───────────────┴────────────────┐
              ↓                                ↓
       perbase_error                    uncalled4_align
       (standard error track)           (signal alignment, 6-16 min)
              ↓                                ↓
       data/perbase_error/         uncalled4_convert_tsv (×2: for + rev)
                                               ↓
                                    data/uncalled4_tsv/
                                               ↓
                                    perbase_signal_deviation
                                               ↓
                                    data/perbase_signal/
                                    (10-col format; col 10 = mean squared deviation
                                     used by phase 3b via Calculate_reactivity.sh -f 10)
```

The per-base signal deviation files have 10 columns. Phase 3b uses `Calculate_reactivity.sh -f 10` to extract the mean squared deviation column (column 10) and compute treatment − control differences. All downstream bedGraph and bigWig steps (phases 4b) receive the resulting 4-column reactivity output and reuse the same scripts (`react_to_bg.sh`, `bg_to_bw.sh`, etc.) with different input/output paths.

### Uncalled4 Known Issues

| Issue | Symptom | Fix |
|-------|---------|-----|
| Non-zero exit on partial DTW failure | Script aborts even though most reads succeeded | `|| true` in scripts + non-empty output check |
| `*` sentinel for DTW failures | TSV column has `*` instead of a number | `na_values=["*"]` in `pd.read_csv` |
| Missing newlines between TSV rows | Row with too many fields | `on_bad_lines="warn"` in `pd.read_csv` |
| pod5/pyarrow deadlock | Hangs indefinitely inside C extension | Pin `pyarrow>=14,<20` (`lib-pod5==0.3.10` does not exist on PyPI; any recent pod5 version is fine) |
| setuptools incompatibility | `pip install uncalled4` build fails | `pip install setuptools==69.5.1` first |
| `-p` processes defaults to 1 | Single-threaded despite Snakemake allocating 8 threads | Always pass `-p "$THREADS"` explicitly |

### dtw.model_diff Sign Convention

`dtw.model_diff` = **model current − observed current** (pA). This is the official uncalled4 definition (predicted minus actual):
- **Positive** value → observed ion current is **lower** than the pore model expects
- **Negative** value → observed current is higher than expected

This sign is consistent across treatment and control samples, so the signal reactivity calculation (treatment_deviation − control_deviation) correctly captures chemical modification-induced changes.

### SLURM Resource Requirements for Uncalled4

From js4007 benchmarks on WGS data:

| Metric | Range |
|--------|-------|
| Wall time per sample | 6–16 minutes |
| Peak RSS | 23–75 GB |
| Effective CPU cores | 6–11 (560–1078% CPU load) |

**Recommended SLURM request:**
- `--mem=80G` (to accommodate worst-case 75 GB RSS)
- `--ntasks=8` or `--cpus-per-task=8` (uncalled4 with `-p 8`)
- `--time=0:30:00` (30 minutes is sufficient for most WGS samples)

The `uncalled4_align` Snakemake rule uses `threads: 8` which reserves 64 GB on Roar's standard partition (8 GB/core), which is adequate for most samples but may be marginal for the worst-case 75 GB peak.
