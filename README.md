# NanoPrint_toolkit

# Overview

NanoPrint_toolkit is a Snakemake-based pipeline for analyzing chemical footprinting data from long-read (Oxford Nanopore) sequencing. The pipeline processes raw sequencing reads through alignment, per-base error quantification, reactivity calculation (treatment minus control), and downstream analyses including significance-filtered bigWig generation, reactive nucleotide density calculations, and feature annotation. The workflow is optimized for parallelization through chromosome-level splitting of data.

However, the toolkit is composed of a series of scripts found in workflow/scripts. Thus, a user can use the pipeline as intended with snakemake, or use individal scripts as documented below.

# Inputs

## Raw Sequencing Reads

### Description
Oxford Nanopore sequencing reads for treatment and control samples. Treatment samples are typically treated with a chemical probe (e.g., permanganate) while control samples are untreated.

### Format
Standard FASTQ (gzipped) or unaligned BAM format.

### Location
`raw_data/{sample}.fastq.gz` or `raw_data/{sample}.bam`

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
Per-base error rates calculated from aligned reads, separated by strand. Error rate represents the probability of a sequencing error at each genomic position.

### Format
Tab-delimited, gzipped, no header:
| Column | Name | Description |
|--------|------|-------------|
| 1 | chrom | Chromosome name |
| 2 | position | Genomic position (1-based) |
| 3 | nucleotide | Reference nucleotide (A, C, G, T) |
| 4 | coverage | Read coverage at position |
| 5 | error | Per-base error probability |

### Example
```
chr19_MATERNAL	1	C	45	0.022222
chr19_MATERNAL	2	C	47	0.021277
chr19_MATERNAL	3	T	48	0.020833
chr19_MATERNAL	4	A	52	0.019231
chr19_MATERNAL	5	A	55	0.018182
chr19_MATERNAL	6	C	58	0.017241
chr19_MATERNAL	7	C	60	0.016667
chr19_MATERNAL	8	C	62	0.016129
chr19_MATERNAL	9	T	65	0.015385
chr19_MATERNAL	10	A	67	0.014925
```

### Location
`data/perbase_error/{genome}/{sample}_{strand}.txt.gz`

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

# Dependencies

All dependencies can be installed via conda:

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
| Languages | python (>=3.8), R |
| R packages | ggplot2, dplyr |

# Instructions

## 1. Organize data

Make directories for inputs

```bash
mkdir -p resources/genomes resources/features raw_data
```
Then add the relevant files in to the correct directory.

## 2. Edit CONFIG file

Your CONFIG file should look something like this:

```
#Genomes, place in resources/genomes/
^g	test_genome.fa

#Features to annotate, place in resources/features/
#^f	bed file	n_windows	window_size
^f	g4Discovery.bed	1000	10

# Window size
^w	1000000
^w	10000

#Significance threshold level, 1 is p <= 0.05,  2 is p <= 0.01,  3 is p <= 0.001,  4 is p <= 0.0001.
^s	2

#Relationships between files
#^    Sample          Treatment                 Control
^r	Hsap_HG002_LCL	Hsap_HG002_LCL_Mn04.bam	Hsap_HG002_LCL_CTRL.bam

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
```
The wildcard variables are assigned designated as:

^g The genome you want to map to
^f Features to annotate with optional parameters: bed_file, n_windows (default: 1000), window_size (default: 10)
^w The window files you want in the windows bed files
^s The significance threshold for identifying reactive nucleotides
^r The relationship between sequencing reads
^t Directories containing temporary files (auto-deleted after use)
^igv-bam Generate strand-split BAMs and indices for IGV visualization (flag, no value)
^igv-bigwig Generate coverage bigWig files for each strand-split BAM (flag, no value)

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

Direct execution with Snakmake

```bash
snakemake --cores 4
```

# Pipeline

```
 ┌─────────────────────────────────────────────┐
 │           PHASE 1: MAPPING & QC             │
 └─────────────────────────────────────────────┘
              ┌─────────────────────────────────────────────┐
              │      Raw reads (fastq, unaligned bam)       │
              └─────────────────────────────────────────────┘
                                     │
                    ┌────────────────┼───────────────────┐
                    │                                    │                                    
                    ▼                                    ▼                                    
        ┌───────────────────┐              ┌───────────────────┐                             
        │  1. read_stats    │              │   2. map_reads    │                             
        │  (Read_stats.sh)  │              │  (Map_reads.sh)   │                             
        └───────────────────┘              └─────────┬─────────┘                             
                                                     │                                       
                                                     ▼                                       
                                          ┌───────────────────────┐                          
                                          │ 3. filter_alignments  │                          
                                          │(Filter_alignments.sh) │                          
                                          └─────────┬─────────────┘                          
                                                    │                                        
                       ┌────────────────────────────┼──────────────────────────┐                         
                       │                            │                          │                         
                       ▼                            │                          ▼                         
        ┌───────────────────┐                       │              ┌───────────────────┐               
        │ 4. alignment_stats│                       │              │   5. histograms   │               
        │(Alignment_stats.sh│                       │              │(Make_histograms.sh│               
        └───────────────────┘                       │              └───────────────────┘               
                                                    │                                         
┌─────────────────────────────────────────────┐     │   
│      PHASE 2: PER-BASE ERROR                │     │     
└─────────────────────────────────────────────┘     │     
                                                    │                                        
                                                    ▼                                   
                                          ┌───────────────────────┐                          
                                          │   6. perbase_error    │
                                          │  (perbase_error.sh)   │
                                          └──────────┬────────────┘
                                                     │
                              ┌──────────────────────┼──────────────────────┐
                              │                                             │
                              ▼                                             ▼
                   ┌─────────────────────────┐                  ┌───────────────────┐
                   │ 7. split_perbase_by_chr │                  │  8. correlation   │
                   │    (Split_by_chr.sh)    │                  │ (Correlation.sh)  │
                   │      [CHECKPOINT]       │                  └───────────────────┘
                   └───────────┬─────────────┘       
                               │                     
                               │  (per chromosome)   
                               │                     
┌───────────────────────┐      ┼──────────────────────────────────┐
│ PHASE 3: REACTIVITY   │      │                                  │
└───────────────────────┘      │                                  │
                               │                                  │
                               ▼                                  │
                     ┌─────────────────────────┐                  │
                     │ 9. calculate_reactivity │                  │
                     │(Calculate_reactivity.sh)│                  │
                     └───────────┬─────────────┘                  │
                                 │                                │
              ┌──────────────────┼───────────────────────────┐    │
              │                  │                           │    │
              │                  │                           │    │
┌─────────┐   │                  │                           │    │
│PHASE 4: │   │                  │                           │    │
│OUTPUT   │   │                  │                           │    │
│FORMATS  │   │                  │                           │    │
└─────────┘   │                  │                           │    │
              │                  │                           │    │
              ▼                  ▼                           │    │
┌─────────────────────────┐   ┌─────────────────────────┐    │    │
│10. reactivity_to_bedgraph│  │  12. bedgraph_to_bigwig │    │    │
│   (react_to_bg.sh)      │   │     (bg_to_bw.sh)       │    │    │
└───────────┬─────────────┘   └───────────┬─────────────┘    │    │
            │                             │ (merge           │    │
            ▼                             ▼  chromosomes)    │    │
┌─────────────────────────┐   ┌─────────────────────────┐    │    │
│ 11. reactivity_density  │   │   13. merge_bigwig      │    │    │
│   (react_dens.sh)       │   │   (Merge_bigwig.sh)     │    │    │
└───────────┬─────────────┘   └─────────────────────────┘    │    │
            │ (merge chromosomes)                            │    │
            ▼                                                │    │
┌─────────────────────────┐                                  │    │
│   14. merge_density     │                                  │    │
│   (Merge_density.sh)    │                                  │    │
└─────────────────────────┘                                  │    │
                                                             │    │
┌──────────────────────────────┐                             │    │
│ PHASE 5: FEATURE ANNOTATION  │                             │    │
└──────────────────────────────┘                             │    │
                                                             │    │
                                                             │    │
                              ┌─────────────────────────┐    │    │
                              │15. split_features_by_chr│    │    │
                              │    (Split_by_chr.sh)    │    │    │
                              │      [CHECKPOINT]       │    │    │
                              └───────────┬─────────────┘    │    │
                                          │                  │    │
                                          ▼                  │    │
                              ┌─────────────────────────┐    │    │
                              │  16. annotate_features  │◄───┘    │
                              │ (annotate_features.sh)  │◄────────┘
                              └───────────┬─────────────┘
                                          │ (merge chromosomes)
                                          ▼
                              ┌─────────────────────────┐
                              │  17. merge_annotations  │
                              │ (Merge_annotations.sh)  │
                              └───────────┬─────────────┘
                                          │
                                          ▼
                              ┌───────────────────────────────┐
                              │ 18. average_annotations       │
                              │(average_feature_annotation.sh)│
                              └───────────────────────────────┘
```

# Commands

---

## 1. read_stats

**Description:** Calculate sequencing statistics from raw reads.

**Script:** `workflow/scripts/Read_stats.sh`

**Inputs:**
- `raw_data/{sample}.fastq.gz`

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
- `raw_data/{sample}.fastq.gz`
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

**Description:** Extract histogram data and generate PDF plots from stats files.

**Script:** `workflow/scripts/Make_histograms.sh`

**Inputs:**
- `data/{alignment}/{genome}/{sample}_stats.txt`

**Outputs:**
- `data/{alignment}/{genome}/{sample}_histograms.txt`
- `data/{alignment}/{genome}/{sample}_histograms.pdf`

**Dependencies:**
- gawk
- R (with ggplot2 for PDF output)

**Documentation:**

```
Usage: Make_histograms.sh -i <input_stats.txt> -o <output_histograms.txt> [-p <output_histograms.pdf>]

Extract histogram data from samtools stats output files.

Required arguments:
    -i    Input stats file (from samtools stats / Alignment_stats.sh)
    -o    Output histogram file (tab-delimited)

Optional arguments:
    -p    Output PDF file with histogram plots
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
    Make_histograms.sh -i sample_stats.txt -o sample_histograms.txt -p sample_histograms.pdf
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
- gawk

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
    Column 5: Per-base error probability

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
    -p    Significance level filter (default: 4)
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

## 15. split_features_by_chr

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

## 16. annotate_features

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

## 17. merge_annotations

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

## 18. average_annotations

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

## 19. igv_split_bam

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

## 20. igv_index_bam

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

## 21. igv_coverage_bigwig

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
