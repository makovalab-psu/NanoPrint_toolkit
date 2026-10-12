#!/usr/bin/env python3
"""
perread_signal.py - Per-read, per-nucleotide signal table: one row for every
reference position of every read of one strand.

The single-molecule view of the data. Nothing is averaged or tested here: each row
is the signal one read gave at one position, with the read's name, so that reads
can be followed along their length and the same read can be found again in the
other backend's table and in the per-read modification calls
(perbase_mod_calls.py --per-read), which share read_id, chr, pos and strand.

Output (tab-delimited, gzipped, with a header):
    read_id
    chr
    pos         1-based reference position
    strand      + or -
    nt          reference base at pos (forward strand)
    kmer        the k-mer the level is modelled by, in read orientation (5'->3' as
                the strand went through the pore); NA if it holds a non-ACGT base
    level       uncalled4: dtw.current, mean normalized current of the samples DTW
                           assigned to the position
                remora:    trimmed-mean normalized signal of the base
    level_sd    SD of those samples
    dwell       number of signal samples at the position
    model_diff  uncalled4 only: level minus the pore model's level for the k-mer

The two backends' levels are in different units and are registered to different
bases of the k-mer, so compare patterns between them, not values.

SIZE. One row per aligned base per read, about 90 bytes before compression: fine
for model oligos, plasmids and single loci, not for a genome (js4022, 4.2 Mb at
~1000x, would be ~9.5e9 rows per sample strand). This is why the workflow only
builds it when CONFIG asks with '^per-read'.

This replaces 'uncalled4 convert' for this purpose: that TSV identifies a read only
by a running number, so it cannot be joined to anything, and convert's own memory
use was never solved on deep data. The uncalled4_convert_tsv rule is still there
if its exact output is wanted.

Usage:
    perread_signal.py --backend uncalled4 -i <uncalled4.bam> -g <genome.fa> \
        -s <for|rev> -o <output.tsv.gz>
    perread_signal.py --backend remora -i <moves.bam> --pod5 <pod5_dir> \
        --levels <9mer_levels.txt> -g <genome.fa> -s <for|rev> -o <output.tsv.gz>
"""

import os
import sys
import argparse

import numpy as np

from mod_calling_common import add_backend_args, open_source, code_to_kmer, open_out


def parse_args():
    p = argparse.ArgumentParser(
        description="Per-read, per-nucleotide signal table",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    add_backend_args(p)
    p.add_argument("-s", "--strand", required=True, choices=["for", "rev"],
                   help="Strand to process")
    p.add_argument("-o", "--output", required=True, help="Output file (.tsv.gz)")
    return p.parse_args()


def fmt(x):
    return "NA" if np.isnan(x) else f"{x:.4f}"


def main():
    args = parse_args()

    print("=== Per-read Signal Table ===")
    print(f"Backend:   {args.backend}")
    print(f"Input BAM: {args.bam}")
    print(f"Strand:    {args.strand}")
    print(f"Reads:     {args.half}")
    print(f"Output:    {args.output}")
    print()

    if not os.path.exists(args.bam):
        sys.exit(f"Error: input BAM not found: {args.bam}")
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    source = open_source(args, strand=args.strand, need_clean=False, want_extra=True)
    k = source.kmer_len
    strand_char = "-" if args.strand == "rev" else "+"
    has_model_diff = args.backend == "uncalled4"
    kmer_names = {}
    n_reads = n_rows = 0

    with open_out(args.output) as out:
        out.write("read_id\tchr\tpos\tstrand\tnt\tkmer\tlevel\tlevel_sd\tdwell"
                  + ("\tmodel_diff" if has_model_diff else "") + "\n")
        for read in source:
            n_reads += 1
            sequence, _ = source.refs.get(read.ref)
            sd, dwell = read.extra["level_sd"], read.extra["dwell"]
            diff = read.extra["model_diff"] if has_model_diff else None
            # Rows in reference order within the read, whichever way it was read.
            for j in np.argsort(read.pos, kind="stable"):
                p0 = int(read.pos[j])
                code = int(read.kmer[j])
                if code not in kmer_names:
                    kmer_names[code] = code_to_kmer(code, k) if code >= 0 else "NA"
                nt = sequence[p0] if 0 <= p0 < len(sequence) else "N"
                out.write(
                    f"{read.read_id}\t{read.ref}\t{p0 + 1}\t{strand_char}\t{nt}"
                    f"\t{kmer_names[code]}\t{read.level[j]:.4f}\t{fmt(sd[j])}"
                    f"\t{'NA' if np.isnan(dwell[j]) else int(dwell[j])}"
                    + (f"\t{fmt(diff[j])}" if has_model_diff else "") + "\n")
                n_rows += 1

    source.report()
    source.close()
    print(f"Reads written: {n_reads}")
    print(f"Rows written:  {n_rows}")
    print("Done.")


if __name__ == "__main__":
    main()
