#!/usr/bin/env python3
"""
perbase_mod_calls.py - Per-read modification calls against a canonical k-mer model,
summarised per reference position.

Steps 2 and 3 of direct modification calling (see mod_calling_common.py). Every
observation of one strand is compared with the model built by kmer_signal_model.py:
z = (level - k-mer mean) / k-mer SD, two-sided Gaussian p-value. Unlike the model,
nothing is excluded for alignment errors here - a modified base is exactly what
makes the basecaller err.

The threshold for calling a base modified is set one of two ways:

  --pval P    (default, P = 0.02 as Rembo)  modified where the p-value is below P.
              Trusts the Gaussian null: on unmodified DNA about P of the
              observations are called, more wherever the level distribution has
              heavier tails than a Gaussian.
  --fpr F --null-hist H [H ...]
              modified where the statistic exceeds the cutoff that a fraction F of
              the HELD-OUT CONTROL exceeds. H are histograms written by --hist-out
              when the control's held-out reads were scored against the same model.
              Empirical: the control comes out at F overall by construction,
              whatever the shape of the null. The cutoff is global - one value for
              all positions and both strands.

Output format (tab-delimited, gzipped; columns 1-5 match perbase_error, so
Calculate_reactivity.sh and the phase 4 scripts can read it):
    1. Chromosome
    2. Position (1-based)
    3. Nucleotide
    4. Coverage (reads scored at this position)
    5. Modified fraction (reads called modified / coverage)
    6. Modified reads
    7. Mean z      - signed; which way the current moved
    8. Mean |z|    - threshold-free effect size

Positions whose k-mer is absent from the model, or has fewer than --min-kmer-obs
observations there, are not scored and do not count toward coverage.

--hist-out writes the distribution of the calling statistic over every scored
observation (|z|, or -log10 of the Fisher-combined p when --fisher-lag > 0), in
bins of 0.01. With --hist-out, -o may be omitted.

--per-read additionally writes one row per read per scored position:
    read_id  chr  pos  strand  level  z  p
(p is the Fisher-combined value when --fisher-lag > 0). Large; meant for model
oligos, plasmids and single loci.

Usage:
    perbase_mod_calls.py --backend uncalled4 -i <uncalled4.bam> -g <genome.fa> \
        -m <model.tsv> -s <for|rev> -o <output.txt.gz> [--half B] [--pval 0.02]
    perbase_mod_calls.py ... --half B --hist-out <control_for.hist.tsv>
    perbase_mod_calls.py ... --fpr 0.02 --null-hist <control_for.hist.tsv> \
        <control_rev.hist.tsv> -o <output.txt.gz>
"""

import os
import sys
import math
import argparse

import numpy as np

from mod_calling_common import (add_backend_args, open_source, load_model, z_scores,
                                two_sided_p, z_threshold, fisher_window,
                                PositionAccumulator, open_out, hist_add, write_hist,
                                cutoff_from_hists, HIST_BINS, SMALLEST_P)


def parse_args():
    p = argparse.ArgumentParser(
        description="Per-read modification calls against a canonical k-mer model",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    add_backend_args(p)
    p.add_argument("-m", "--model", required=True,
                   help="K-mer model from kmer_signal_model.py (same backend)")
    p.add_argument("-s", "--strand", required=True, choices=["for", "rev"],
                   help="Strand to process")
    p.add_argument("-o", "--output", help="Output file (.txt.gz)")
    p.add_argument("--pval", type=float,
                   help="A base is modified below this p-value (default: 0.02, as "
                        "Rembo, unless --fpr is given)")
    p.add_argument("--fpr", type=float,
                   help="Instead of --pval: call at the cutoff this fraction of the "
                        "held-out control exceeds. Needs --null-hist.")
    p.add_argument("--null-hist", nargs="+",
                   help="Histograms from --hist-out on the held-out control (same "
                        "model, same --fisher-lag); all are pooled")
    p.add_argument("--hist-out", help="Write the calling statistic's histogram here")
    p.add_argument("--min-kmer-obs", type=int, default=30,
                   help="Do not score k-mers the model saw fewer times than this "
                        "(default: 30)")
    p.add_argument("--fisher-lag", type=int, default=0,
                   help="Combine each p-value with its +/- N neighbours along the read "
                        "by Fisher's method before calling (default: 0 = off; Tombo's "
                        "default is 1)")
    p.add_argument("--per-read", help="Also write the per-read table here (.tsv.gz)")
    p.add_argument("-c", "--min-cov", type=int, default=1,
                   help="Minimum coverage to emit a position (default: 1)")
    args = p.parse_args()
    if args.pval is not None and args.fpr is not None:
        p.error("--pval and --fpr are alternatives; give one")
    if args.fpr is not None and not args.null_hist:
        p.error("--fpr needs --null-hist (histograms of the held-out control)")
    if args.null_hist and args.fpr is None:
        p.error("--null-hist is only used with --fpr")
    if not args.output and not args.hist_out:
        p.error("nothing to write: give -o and/or --hist-out")
    if args.fpr is None and args.pval is None:
        args.pval = 0.02
    return args


def main():
    args = parse_args()

    # The statistic a call is made on: |z|, or -log10 p once p-values have been
    # combined along the read. Both grow with evidence, so a call is stat > cutoff.
    statistic = "neglog10p" if args.fisher_lag else "absz"
    hist_meta = {"statistic": statistic, "fisher_lag": args.fisher_lag,
                 "model": args.model}
    if args.fpr is not None:
        cutoff, n_null = cutoff_from_hists(args.null_hist, args.fpr, hist_meta)
        how = (f"false-positive rate {args.fpr} of {n_null} held-out control "
               f"observations")
    elif args.fisher_lag:
        cutoff = -math.log10(args.pval)
        how = f"p < {args.pval}"
    else:
        cutoff = z_threshold(args.pval)
        how = f"p < {args.pval}"

    print("=== Per-base Modification Calls ===")
    print(f"Backend:   {args.backend}")
    print(f"Input BAM: {args.bam}")
    print(f"Model:     {args.model}")
    print(f"Strand:    {args.strand}")
    print(f"Reads:     {args.half}")
    print(f"Threshold: {how}  ->  {statistic} > {cutoff:.4f}"
          + (f", Fisher lag {args.fisher_lag}" if args.fisher_lag else ""))
    print(f"Output:    {args.output or '(histogram only)'}")
    print()

    for path in (args.bam, args.model):
        if not os.path.exists(path):
            sys.exit(f"Error: input not found: {path}")
    for path in (args.output, args.hist_out, args.per_read):
        if path:
            os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)

    meta, m_kmer, m_n, m_mean, m_sd = load_model(args.model)
    source = open_source(args, strand=args.strand, need_clean=False)

    # Levels from the two backends differ in units and in which base of the k-mer
    # they belong to; a model from the other backend would score without complaint
    # and mean nothing.
    context = f"{source.context[0]},{source.context[1]}"
    for key, have in (("backend", args.backend), ("kmer_len", str(source.kmer_len)),
                      ("kmer_context", context), ("levels", args.levels or "none")):
        if meta.get(key) != have:
            sys.exit(f"Error: model {args.model} was built with {key}={meta.get(key)}, "
                     f"but this run has {key}={have}.")

    need_p = bool(args.fisher_lag or args.per_read)
    strand_char = "-" if args.strand == "rev" else "+"
    written = 0
    n_obs = n_scored = n_called = 0
    hist = np.zeros(HIST_BINS + 1, dtype=np.int64)

    out = open_out(args.output) if args.output else None

    def emit(ref, positions, cov, n_mod, sum_z, sum_abs):
        nonlocal written
        if out is None:
            return
        sequence, _ = source.refs.get(ref)
        for k in range(positions.size):
            c = int(cov[k])
            if c < args.min_cov:
                continue
            p0 = int(positions[k])
            nt = sequence[p0] if 0 <= p0 < len(sequence) else "N"
            if nt not in "ACGTN":
                nt = "N"
            out.write(f"{ref}\t{p0 + 1}\t{nt}\t{c}\t{n_mod[k] / c:.6f}"
                      f"\t{int(n_mod[k])}\t{sum_z[k] / c:.6f}"
                      f"\t{sum_abs[k] / c:.6f}\n")
            written += 1

    acc = PositionAccumulator(emit)
    per_read = open_out(args.per_read) if args.per_read else None
    if per_read:
        per_read.write("read_id\tchr\tpos\tstrand\tlevel\tz\tp\n")

    for read in source:
        n_obs += read.pos.size
        z = z_scores(read.kmer, read.level, m_kmer, m_n, m_mean, m_sd,
                     args.min_kmer_obs)
        ok = ~np.isnan(z)
        pos, z, level = read.pos[ok], z[ok], read.level[ok]
        if need_p:
            p = two_sided_p(z)
        if args.fisher_lag:
            p = fisher_window(pos, p, args.fisher_lag)
            stat = -np.log10(np.maximum(p, SMALLEST_P))
        else:
            # |z| past its threshold is the same test as p below its threshold,
            # without evaluating erfc for every observation.
            stat = np.abs(z)
        called = stat > cutoff
        if args.hist_out:
            hist_add(hist, stat)
        n_scored += pos.size
        n_called += int(called.sum())
        acc.set_ref(read.ref)
        acc.add(pos, z, called, read.start)
        if per_read:
            for j in range(pos.size):
                per_read.write(f"{read.read_id}\t{read.ref}\t{int(pos[j]) + 1}"
                               f"\t{strand_char}\t{level[j]:.4f}\t{z[j]:.4f}"
                               f"\t{p[j]:.4g}\n")
    acc.close()
    if per_read:
        per_read.close()
    if out:
        out.close()
    if args.hist_out:
        write_hist(args.hist_out, hist, hist_meta)

    source.report()
    source.close()

    print(f"Observations with a level: {n_obs}")
    print(f"Scored against the model:  {n_scored}")
    if n_scored:
        print(f"Called modified:           {n_called} ({n_called / n_scored:.4%})")
    if out:
        print(f"Positions written:         {written}")
    if n_obs > 0 and n_scored == 0:
        sys.exit("Error: observations were read but none could be scored - no k-mer "
                 "in this sample is in the model with enough observations. Check that "
                 "the model was built on the same reference, and --min-kmer-obs.")
    print("Done.")


if __name__ == "__main__":
    main()
