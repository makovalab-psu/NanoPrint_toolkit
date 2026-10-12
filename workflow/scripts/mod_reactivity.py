#!/usr/bin/env python3
"""
mod_reactivity.py - Modification-rate reactivity: treatment minus control modified
fraction, per position.

Both inputs are perbase_mod_calls.py outputs scored against the SAME k-mer model.
The control's modified fraction is the background the test produces on unmodified
DNA - about the p-value threshold by construction, more where the level
distribution is not Gaussian - so the difference is what the treatment added.

Unlike Calculate_reactivity.sh this takes whole-genome files (every contig in one
file) and merges them by contig and position, so no per-chromosome split is needed.
Both inputs must list contigs in the order of the genome .fai, which they do: they
are written in BAM coordinate order.

Output format (tab-delimited, gzipped; columns 1-4 match the reactivity files):
    1. Chromosome
    2. Position (1-based)
    3. Nucleotide
    4. Reactivity (treatment - control modified fraction)
    5. Treatment coverage
    6. Control coverage
    7. Treatment modified fraction
    8. Control modified fraction

Only positions scored in both samples at >= --min-cov are written.

Usage:
    mod_reactivity.py -p <treatment.txt.gz> -m <control.txt.gz> -g <genome.fa.fai> \
        -o <output.txt.gz> [-c 10]
"""

import sys
import gzip
import argparse


def parse_args():
    p = argparse.ArgumentParser(
        description="Treatment minus control modified fraction per position",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    p.add_argument("-p", "--treatment", required=True)
    p.add_argument("-m", "--control", required=True)
    p.add_argument("-g", "--fai", required=True, help="Genome index (.fai): contig order")
    p.add_argument("-o", "--output", required=True)
    p.add_argument("-c", "--min-cov", type=int, default=10,
                   help="Minimum coverage in BOTH samples (default: 10)")
    return p.parse_args()


def rows(path, rank):
    """Yield ((contig rank, position), fields), checking the file is in that order."""
    last = None
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] not in rank:
                sys.exit(f"Error: contig {f[0]} in {path} is not in the genome index.")
            key = (rank[f[0]], int(f[1]))
            if last is not None and key <= last:
                sys.exit(f"Error: {path} is not sorted by contig and position at "
                         f"{f[0]}:{f[1]}; the merge would silently skip rows.")
            last = key
            yield key, f


def main():
    args = parse_args()
    with open(args.fai) as fh:
        rank = {line.split("\t", 1)[0]: i for i, line in enumerate(fh)}

    treatment = rows(args.treatment, rank)
    control = rows(args.control, rank)
    t = next(treatment, None)
    c = next(control, None)
    n_t = n_c = written = 0
    with gzip.open(args.output, "wt") as out:
        while t is not None and c is not None:
            if t[0] < c[0]:
                n_t += 1
                t = next(treatment, None)
            elif c[0] < t[0]:
                n_c += 1
                c = next(control, None)
            else:
                ft, fc = t[1], c[1]
                if int(ft[3]) >= args.min_cov and int(fc[3]) >= args.min_cov:
                    out.write(f"{ft[0]}\t{ft[1]}\t{ft[2]}"
                              f"\t{float(ft[4]) - float(fc[4]):.6f}"
                              f"\t{ft[3]}\t{fc[3]}\t{ft[4]}\t{fc[4]}\n")
                    written += 1
                n_t += 1
                n_c += 1
                t = next(treatment, None)
                c = next(control, None)
        n_t += sum(1 for _ in treatment) + (t is not None)
        n_c += sum(1 for _ in control) + (c is not None)

    print(f"Treatment positions: {n_t}")
    print(f"Control positions:   {n_c}")
    print(f"Positions written:   {written} (both >= {args.min_cov}x)")
    if n_t and n_c and written == 0:
        sys.stderr.write("Warning: both inputs have positions but none passed in both. "
                         "Check --min-cov against the per-strand coverage.\n")


if __name__ == "__main__":
    main()
