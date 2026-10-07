# Phase 6: IGV Export
# Rules for generating strand-split BAMs, indices, and coverage bigWigs


def get_igv_source_bam(wildcards):
    """Get the source BAM file for an IGV export rule.

    A ^u sample has no BAM under data/; the supplied Uncalled4 BAM is the
    filtered alignment, and rule all asks only for its filtered_alignments
    export (see igv_sources()).
    """
    if wildcards.igv_source == "filtered_alignments" and has_uncalled4(wildcards.raw_sample):
        return get_uncalled4_bam(wildcards.raw_sample)
    return f"{IGV_SOURCE_DIRS[wildcards.igv_source]}/{wildcards.genome}/{wildcards.raw_sample}.bam"


rule igv_split_bam:
    """Split BAM into forward and reverse strand reads for IGV visualization."""
    input:
        bam=get_igv_source_bam
    output:
        bam="results/igv/{igv_source}/{genome}/{raw_sample}_{strand}.bam"
    params:
        flag=lambda wc: "-F 0x10" if wc.strand == "for" else "-f 0x10"
    log:
        "logs/igv/split_bam/{igv_source}/{genome}/{raw_sample}_{strand}.log"
    benchmark:
        "benchmarks/phase6/igv_split_bam/{igv_source}/{genome}/{raw_sample}_{strand}.tsv"
    wildcard_constraints:
        strand="for|rev",
        igv_source="all_alignments|filtered_alignments"
    shell:
        """
        samtools view -b -h {params.flag} {input.bam} > {output.bam} 2> {log}
        """


rule igv_index_bam:
    """Create BAM index for IGV visualization."""
    input:
        bam="results/igv/{igv_source}/{genome}/{raw_sample}_{strand}.bam"
    output:
        bai="results/igv/{igv_source}/{genome}/{raw_sample}_{strand}.bam.bai"
    log:
        "logs/igv/index_bam/{igv_source}/{genome}/{raw_sample}_{strand}.log"
    benchmark:
        "benchmarks/phase6/igv_index_bam/{igv_source}/{genome}/{raw_sample}_{strand}.tsv"
    wildcard_constraints:
        strand="for|rev",
        igv_source="all_alignments|filtered_alignments"
    shell:
        """
        samtools index {input.bam} 2> {log}
        """


rule igv_coverage_bigwig:
    """Generate coverage bigWig from strand-split BAM for IGV visualization."""
    input:
        bam="results/igv/{igv_source}/{genome}/{raw_sample}_{strand}.bam",
        bai="results/igv/{igv_source}/{genome}/{raw_sample}_{strand}.bam.bai",
        fai="resources/genomes/{genome}.fa.fai"
    output:
        bw="results/igv/{igv_source}/{genome}/{raw_sample}_{strand}.bw"
    log:
        "logs/igv/coverage_bigwig/{igv_source}/{genome}/{raw_sample}_{strand}.log"
    benchmark:
        "benchmarks/phase6/igv_coverage_bigwig/{igv_source}/{genome}/{raw_sample}_{strand}.tsv"
    wildcard_constraints:
        strand="for|rev",
        igv_source="all_alignments|filtered_alignments"
    shell:
        """
        TMP_DIR=$(mktemp -d)
        trap "rm -rf $TMP_DIR" EXIT

        # Count from the INDEX. Do NOT revert this to `samtools view | head -1`:
        # under the `set -euo pipefail` snakemake applies, head closes the pipe,
        # samtools dies of SIGPIPE, and every large BAM is misreported as empty.
        # See claude.md, "samtools view | head -1 is a silent lie under pipefail".
        N_READS=$(samtools idxstats {input.bam} | awk '{{n += $3 + $4}} END {{print n + 0}}')
        if [[ "$N_READS" -eq 0 ]]; then
            echo "No reads in {input.bam} — creating empty bigWig" | tee {log}
            touch {output.bw}
        else
            echo "Reads in {input.bam}: $N_READS" | tee {log}
            bedtools genomecov -ibam {input.bam} -bg | \
                sort -k1,1 -k2,2n > "$TMP_DIR/coverage.bg"
            cut -f1,2 {input.fai} > "$TMP_DIR/chrom.sizes"
            # tee -a, not tee: a second plain `tee {log}` truncates the file and
            # throws away the read count logged above.
            bedGraphToBigWig "$TMP_DIR/coverage.bg" "$TMP_DIR/chrom.sizes" {output.bw} \
                2>&1 | tee -a {log}
        fi
        """
