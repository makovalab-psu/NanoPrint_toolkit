#!/bin/bash

# Demux_bam.sh - Split a barcoded BAM into per-sample BAMs using the dorado BC:Z: tag
#
# Input is typically the Uncalled4 BAM from `nanoprint preprocess -k <kit>`: every read
# carries BC:Z:<kit>_barcodeNN (dorado --kit-name), except unclassified reads, which
# carry no BC tag. A sample sheet maps barcodes to sample names; several barcodes may
# share a name, which pools them into one BAM. Reads with no BC tag, or a barcode not in
# the sheet, go to unclassified.bam.
#
# The split is a single streaming pass (samtools view | awk) with one `samtools view -b`
# pipe per sample, all opened up front so every sample BAM exists even if it gets 0
# reads. samtools, not `dorado demux`, writes the records: samtools keeps every tag,
# including the Uncalled4 DTW tags and MM/ML, while whether dorado demux does is
# unverified. A coordinate-sorted input gives coordinate-sorted outputs, so they are
# indexed directly with no re-sort.
#
# The same pass counts reads and bases per barcode x run (run = read group ID up to
# its first underscore, which is the MinKNOW run ID), so pooled barcodes and runs stay
# separable in demux_counts.tsv.
#
# Usage: Demux_bam.sh -i <input.bam> -s <samplesheet.tsv> -o <out_dir> [-t <threads>]

set -euo pipefail

usage() {
    cat << EOF
Usage: $(basename "$0") -i <input.bam> -s <samplesheet.tsv> -o <out_dir> [-t <threads>]

Split a BAM into per-sample BAMs by the dorado barcode tag (BC:Z:).

Required arguments:
    -i    Input BAM (coordinate-sorted, e.g. from nanoprint preprocess -k <kit>)
    -s    Sample sheet: tab-separated, two columns, barcode<TAB>sample.
            barcode  barcodeNN or <kit>_barcodeNN (e.g. barcode01, SQK-RBK114-24_barcode01)
            sample   output name; letters, digits, '.', '_' and '-' only
          Several barcodes may share a sample name (they are pooled). Blank lines and
          lines starting with '#' are ignored. 'unclassified' is reserved.
    -o    Output directory. Must not already contain .bam files.

Optional arguments:
    -t    Threads for BAM decoding (default: 4)
    -h    Show this help message

Outputs (in <out_dir>):
    <sample>.bam + .bai       one per sample name in the sheet
    unclassified.bam + .bai   reads with no BC tag or a barcode not in the sheet
    demux_counts.tsv          sample, barcode, run_id, reads, bases
    samplesheet.normalized.tsv  the sheet as applied (barcodeNN<TAB>sample)
EOF
    exit 1
}

INPUT=""
SHEET=""
OUT_DIR=""
THREADS=4

while getopts "i:s:o:t:h" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        s) SHEET="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ -z "$INPUT" || -z "$SHEET" || -z "$OUT_DIR" ]]; then
    echo "Error: Missing required arguments" >&2
    usage
fi
[[ -f "$INPUT" ]] || { echo "Error: Input BAM not found: $INPUT" >&2; exit 1; }
[[ -f "$SHEET" ]] || { echo "Error: Sample sheet not found: $SHEET" >&2; exit 1; }
command -v samtools &>/dev/null || { echo "Error: samtools not found in PATH" >&2; exit 1; }

mkdir -p "$OUT_DIR"
if compgen -G "${OUT_DIR}/*.bam" > /dev/null; then
    echo "Error: $OUT_DIR already contains .bam files. Remove them or choose another -o." >&2
    exit 1
fi

TMP_DIR="${OUT_DIR}/tmp_demux_$$"
mkdir -p "$TMP_DIR"
trap "rm -rf '${TMP_DIR}'" EXIT

# ---------------------------------------------------------------------------
# Validate and normalize the sample sheet -> barcodeNN<TAB>sample
# ---------------------------------------------------------------------------
NORM_SHEET="${OUT_DIR}/samplesheet.normalized.tsv"
awk -F'\t' -v OFS='\t' '
    /^[ \t\r]*(#|$)/ { next }
    {
        sub(/\r$/, "")
        if (NF != 2) { printf "Error: sample sheet line %d: expected 2 tab-separated columns, got %d: %s\n", NR, NF, $0 > "/dev/stderr"; bad = 1; next }
        bc = $1; sample = $2
        if (!match(bc, /barcode[0-9]+$/)) { printf "Error: sample sheet line %d: barcode \"%s\" does not end in barcodeNN\n", NR, bc > "/dev/stderr"; bad = 1; next }
        bc = substr(bc, RSTART)
        if (sample !~ /^[A-Za-z0-9._-]+$/) { printf "Error: sample sheet line %d: sample name \"%s\" has characters other than letters, digits, . _ -\n", NR, sample > "/dev/stderr"; bad = 1; next }
        if (sample == "unclassified") { printf "Error: sample sheet line %d: \"unclassified\" is reserved\n", NR > "/dev/stderr"; bad = 1; next }
        if ((bc in seen) && seen[bc] != sample) { printf "Error: sample sheet line %d: %s is assigned to both %s and %s\n", NR, bc, seen[bc], sample > "/dev/stderr"; bad = 1; next }
        if (!(bc in seen)) { seen[bc] = sample; print bc, sample }
    }
    END { exit bad }
' "$SHEET" > "$NORM_SHEET" || { echo "Error: invalid sample sheet: $SHEET" >&2; exit 1; }

[[ -s "$NORM_SHEET" ]] || { echo "Error: sample sheet has no barcode rows: $SHEET" >&2; exit 1; }

HEADER_SAM="${TMP_DIR}/header.sam"
samtools view -H "$INPUT" > "$HEADER_SAM"
SORTED=false
grep -q $'^@HD\t.*SO:coordinate' "$HEADER_SAM" && SORTED=true

echo "=== Demultiplex by BC tag ==="
echo "Input:        $INPUT"
echo "Sample sheet: $SHEET"
echo "Output dir:   $OUT_DIR"
echo "Coordinate-sorted input: $SORTED"
echo ""
echo "Barcode -> sample:"
sed 's/^/    /' "$NORM_SHEET"
echo ""

# ---------------------------------------------------------------------------
# Single pass: route each record to its sample's BAM, count per barcode x run
# ---------------------------------------------------------------------------
COUNTS_RAW="${TMP_DIR}/counts_raw.tsv"
samtools view -@ "$THREADS" "$INPUT" \
  | awk -F'\t' -v OFS='\t' -v sheet="$NORM_SHEET" -v hdr="$HEADER_SAM" \
        -v dir="$OUT_DIR" -v counts="$COUNTS_RAW" '
    BEGIN {
        while ((getline line < hdr) > 0) header = header line "\n"
        close(hdr)
        n = 0
        while ((getline line < sheet) > 0) {
            split(line, f, "\t")
            sample_of[f[1]] = f[2]
            if (!(f[2] in cmd)) names[++n] = f[2]
            cmd[f[2]] = "samtools view -b -@ 2 -o \"" dir "/" f[2] ".bam\" -"
        }
        close(sheet)
        names[++n] = "unclassified"
        cmd["unclassified"] = "samtools view -b -@ 2 -o \"" dir "/unclassified.bam\" -"
        for (i = 1; i <= n; i++) printf "%s", header | cmd[names[i]]
    }
    {
        bc = "unclassified"; run = "NA"
        for (i = 12; i <= NF; i++) {
            tag = substr($i, 1, 5)
            if (tag == "BC:Z:") {
                if (match($i, /barcode[0-9]+$/)) bc = substr($i, RSTART)
            } else if (tag == "RG:Z:") {
                run = substr($i, 6); sub(/_.*/, "", run)
            }
        }
        sample = (bc in sample_of) ? sample_of[bc] : "unclassified"
        print | cmd[sample]
        key = sample SUBSEP bc SUBSEP run
        reads[key]++
        if ($10 != "*") bases[key] += length($10)
    }
    END {
        for (i = 1; i <= n; i++) close(cmd[names[i]])
        for (key in reads) {
            split(key, k, SUBSEP)
            printf "%s\t%s\t%s\t%d\t%.0f\n", k[1], k[2], k[3], reads[key], bases[key] > counts
        }
        close(counts)
    }'

COUNTS="${OUT_DIR}/demux_counts.tsv"
{
    printf "sample\tbarcode\trun_id\treads\tbases\n"
    if [[ -s "$COUNTS_RAW" ]]; then sort -t$'\t' -k1,1 -k2,2 -k3,3 "$COUNTS_RAW"; fi
} > "$COUNTS"

# ---------------------------------------------------------------------------
# Sort (only if the input was not coordinate-sorted) and index
# ---------------------------------------------------------------------------
OUT_BAMS=()
while IFS= read -r sample; do
    OUT_BAMS+=("${OUT_DIR}/${sample}.bam")
done < <(cut -f2 "$NORM_SHEET" | awk '!seen[$0]++'; echo "unclassified")

for bam in "${OUT_BAMS[@]}"; do
    [[ -f "$bam" ]] || { echo "Error: expected output was not written: $bam" >&2; exit 1; }
    if [[ "$SORTED" != true ]]; then
        samtools sort -@ "$THREADS" -T "${TMP_DIR}/sort_tmp" -o "${TMP_DIR}/sorted.bam" "$bam"
        mv "${TMP_DIR}/sorted.bam" "$bam"
    fi
    samtools index "$bam"
done

# ---------------------------------------------------------------------------
# Check nothing was lost: output records must add up to input records
# ---------------------------------------------------------------------------
count_records() {
    samtools idxstats "$1" 2>/dev/null | awk '{ n += $3 + $4 } END { print n + 0 }'
}
if [[ -f "${INPUT}.bai" || -f "${INPUT%.bam}.bai" || -f "${INPUT}.csi" ]]; then
    N_IN=$(count_records "$INPUT")
else
    N_IN=$(samtools view -c -@ "$THREADS" "$INPUT")
fi
N_OUT=0
echo "Records per output BAM:"
for bam in "${OUT_BAMS[@]}"; do
    n=$(count_records "$bam")
    N_OUT=$((N_OUT + n))
    printf "    %-45s %s\n" "$(basename "$bam")" "$n"
    [[ "$n" -eq 0 && "$(basename "$bam")" != unclassified.bam ]] \
        && echo "    Warning: $(basename "$bam") received 0 reads — check its barcodes in the sample sheet" >&2
done
N_COUNTED=$(awk -F'\t' 'NR > 1 { n += $4 } END { print n + 0 }' "$COUNTS")

echo ""
echo "Input records:   $N_IN"
echo "Output records:  $N_OUT"
echo "Counted records: $N_COUNTED"
if [[ "$N_OUT" -ne "$N_IN" || "$N_COUNTED" -ne "$N_IN" ]]; then
    echo "Error: record counts do not match — outputs are incomplete. Do not use them." >&2
    exit 1
fi

echo ""
echo "Done. Counts per barcode x run: $COUNTS"
