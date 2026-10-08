#!/usr/bin/env python3
"""
kmer_signal_model.py - Canonical per-k-mer signal model from an untreated sample.

Step 1 of direct modification calling (see mod_calling_common.py for the method and
its sources). Pools the signal level of every k-mer observation that aligned with no
mismatch or indel, both strands together, removes observations more than --mad
median absolute deviations from the k-mer median, and writes n / mean / SD per k-mer.

Output (tab-delimited, '#key=value' metadata lines, then a header):
    kmer  n  n_outlier  mean  sd  median  mad
n is what mean and sd were computed from; n_outlier is what the MAD filter removed;
median and mad are of the observations before filtering. K-mers are in read
orientation. Only k-mers that were observed are listed - on a short reference that
is a few hundred of the 4^k, and each occurs at about one position, so the "k-mer"
model is in effect a per-position, per-strand model.

Several BAMs (-i a.bam b.bam ...) are pooled into one model, for a control pooled
across replicates or experiments; the remora backend then takes one --pod5 per BAM.
With --max-obs-per-kmer the pooled sample is uniform over observations, so a deeper
BAM contributes proportionally more.

Usage:
    kmer_signal_model.py --backend uncalled4 -i <uncalled4.bam> -g <genome.fa> \
        -o <model.tsv> [--half A]
    kmer_signal_model.py --backend remora -i <moves.bam> --pod5 <pod5_dir> \
        -g <genome.fa> -o <model.tsv> [--levels 9mer_levels_v1.txt] [--half A]
"""

import os
import sys
import argparse

from mod_calling_common import (add_backend_args, open_source, KmerModelBuilder,
                                write_model)


def parse_args():
    p = argparse.ArgumentParser(
        description="Canonical per-k-mer signal model from an untreated sample",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    add_backend_args(p, multi=True)
    p.add_argument("-o", "--output", required=True, help="Output model (.tsv)")
    p.add_argument("--mad", type=float, default=15.0,
                   help="Drop observations more than this many median absolute "
                        "deviations from the k-mer median (default: 15, as Rembo)")
    p.add_argument("--max-obs-per-kmer", type=int, default=20000,
                   help="Uniform random sample kept per k-mer; bounds memory on deep "
                        "data (default: 20000; 0 keeps everything)")
    p.add_argument("--seed", type=int, default=1, help="Seed for that sample (default: 1)")
    p.add_argument("--clean-flank", type=int,
                   help="Use an observation only if the read matches the reference "
                        "within this many bases either side of it (default: k-1, i.e. "
                        "every k-mer that can contain the position). Smaller keeps more "
                        "observations; -1 turns the error mask off and uses all of them.")
    return p.parse_args()


def main():
    args = parse_args()

    print("=== K-mer Signal Model ===")
    print(f"Backend:   {args.backend}")
    print(f"Input BAM: {' '.join(args.bam)}")
    print(f"Genome:    {args.genome}")
    print(f"Reads:     {args.half}")
    print(f"Output:    {args.output}")
    print()

    for bam in args.bam:
        if not os.path.exists(bam):
            sys.exit(f"Error: input BAM not found: {bam}")
    pod5s = args.pod5 or [None] * len(args.bam)
    if len(pod5s) != len(args.bam):
        sys.exit(f"Error: {len(args.bam)} BAMs but {len(pod5s)} --pod5 paths; the "
                 f"remora backend needs one per BAM, in the same order.")
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    builder = KmerModelBuilder(max_obs=args.max_obs_per_kmer, seed=args.seed)
    n_obs = n_candidates = 0
    for bam, pod5 in zip(args.bam, pod5s):
        sys.stderr.write(f"--- {bam}\n")
        source = open_source(args, strand=None, need_clean=True, bam=bam, pod5=pod5)
        for read in source:
            n_obs += read.pos.size
            use = read.clean & (read.kmer >= 0)
            builder.add(read.kmer[use], read.level[use])
        # Report before closing: the read counts are the first thing needed if
        # anything about the source goes wrong, including in close() itself.
        source.report()
        n_candidates += source.counts["candidate reads"]
        source.close()

    table = builder.finalize(args.mad)
    meta = {
        "backend": args.backend,
        "kmer_len": source.kmer_len,
        "kmer_context": f"{source.context[0]},{source.context[1]}",
        "bam": ",".join(args.bam),
        "half": args.half,
        "mad": args.mad,
        "max_obs_per_kmer": args.max_obs_per_kmer,
        "clean_flank": "default" if args.clean_flank is None else args.clean_flank,
        "levels": args.levels or "none",
    }
    write_model(args.output, meta, table)

    print(f"Observations with a level:           {n_obs}")
    print(f"Observations in error-free k-mers:   {builder.seen}")
    print(f"K-mers modelled:                     {table['kmer'].size} "
          f"(k = {source.kmer_len})")
    if table["kmer"].size:
        print(f"Removed by the {args.mad:g}-MAD filter:        "
              f"{int(table['n_outlier'].sum())}")
    if builder.seen == 0 and n_candidates == 0:
        # No reads at all: the sample does not contain this reference (a B-DNA
        # library on the PolyT reverse complement). Same convention as the phase 0
        # scripts - a warning and an empty result, not a failed workflow.
        print("Warning: no reads on this reference; the model is empty. Samples "
              "scored against it will come out empty.")
    elif builder.seen == 0:
        # Reads were decoded yet none had an error-free k-mer. That is not an empty
        # sample, it is something wrong, and an empty model would hide it.
        sys.exit(f"Error: {n_candidates} reads were read but none had an error-free "
                 f"k-mer observation; the model is empty.")
    print("Done.")


if __name__ == "__main__":
    main()
