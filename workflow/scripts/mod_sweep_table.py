#!/usr/bin/env python3
"""
mod_sweep_table.py - Tabulate a modification-calling parameter sweep.

Reads the manifest written by workflow/sweeps/mod_sweep.smk (one row per parameter
set x genome x relationship x control model x strand, naming the treated and the
control perbase_mod_calls.py outputs) and writes two tables.

positions (gzipped TSV) - one row per manifest row x reference position:
    the manifest's parameter and sample columns, then
    chr, pos, nt,
    cov_treated, cov_control          reads scored
    frac_treated, frac_control        fraction of reads called modified
    reactivity                        frac_treated - frac_control
    mean_z_treated, mean_z_control    signed mean z
    mean_absz_treated, mean_absz_control
  A position present in only one sample is kept, with the other side empty.

summary (TSV) - one row per manifest row x contig, over positions where both samples
have at least --min-cov reads:
    n_positions
    cov_treated_median, cov_control_median
    background            mean frac_control - the false-positive rate actually
                          realised on the held-out control
    frac_treated_mean
    reactivity_mean, reactivity_max, pos_of_max, nt_of_max
    reactivity_mean_T, reactivity_mean_nonT, n_T
                          split by the reference base. Permanganate reacts mainly
                          with T, but one modified T shifts the current at its
                          neighbours too, so non-T is not a clean negative.
    reactivity_sd         spread of reactivity across positions
    n_reactive            positions with reactivity > 3 x background

Nothing here chooses between parameter sets; it lays them out to be compared.

Usage:
    mod_sweep_table.py -m manifest.tsv -p positions.tsv.gz -s summary.tsv [-c 10]
"""

import os
import sys
import csv
import gzip
import argparse
import statistics


def parse_args():
    p = argparse.ArgumentParser(
        description="Tabulate a modification-calling parameter sweep",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    p.add_argument("-m", "--manifest", required=True)
    p.add_argument("-p", "--positions", required=True, help="Output (.tsv.gz)")
    p.add_argument("-s", "--summary", required=True, help="Output (.tsv)")
    p.add_argument("-c", "--min-cov", type=int, default=10,
                   help="Coverage both samples need for a position to enter the "
                        "summary (default: 10)")
    return p.parse_args()


def read_calls(path):
    """perbase_mod_calls.py output -> {(chr, pos): (nt, cov, frac, mean_z, mean_absz)},
    plus the contig order of the file."""
    calls, order = {}, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] not in order:
                order.append(f[0])
            calls[(f[0], int(f[1]))] = (f[2], int(f[3]), float(f[4]), float(f[6]),
                                        float(f[7]))
    return calls, order


def fmt(x):
    return "" if x is None else (f"{x:.6f}" if isinstance(x, float) else str(x))


def main():
    args = parse_args()
    with open(args.manifest) as fh:
        manifest = list(csv.DictReader(fh, delimiter="\t"))
    if not manifest:
        sys.exit(f"Error: manifest {args.manifest} has no rows.")
    meta_cols = [c for c in manifest[0] if c not in ("treatment_file", "control_file")]

    pos_cols = ["chr", "pos", "nt", "cov_treated", "cov_control", "frac_treated",
                "frac_control", "reactivity", "mean_z_treated", "mean_z_control",
                "mean_absz_treated", "mean_absz_control"]
    sum_cols = ["chr", "n_positions", "cov_treated_median", "cov_control_median",
                "background", "frac_treated_mean", "reactivity_mean", "reactivity_max",
                "pos_of_max", "nt_of_max", "reactivity_mean_T", "reactivity_mean_nonT",
                "n_T", "reactivity_sd", "n_reactive"]

    cache = {}

    def calls_for(path):
        if path not in cache:
            if not os.path.exists(path):
                sys.exit(f"Error: {path} is listed in the manifest but does not exist.")
            cache[path] = read_calls(path)
        return cache[path]

    n_rows = n_empty = 0
    with gzip.open(args.positions, "wt") as pos_out, open(args.summary, "w") as sum_out:
        pos_out.write("\t".join(meta_cols + pos_cols) + "\n")
        sum_out.write("\t".join(meta_cols + sum_cols) + "\n")
        for unit in manifest:
            meta = [unit[c] for c in meta_cols]
            treated, t_order = calls_for(unit["treatment_file"])
            control, c_order = calls_for(unit["control_file"])
            contigs = t_order + [c for c in c_order if c not in t_order]
            if not treated and not control:
                n_empty += 1
            per_contig = {c: [] for c in contigs}
            for key in sorted(set(treated) | set(control),
                              key=lambda k: (contigs.index(k[0]), k[1])):
                t, c = treated.get(key), control.get(key)
                nt = (t or c)[0]
                react = t[2] - c[2] if t and c else None
                pos_out.write("\t".join(meta + [
                    key[0], str(key[1]), nt,
                    fmt(t[1] if t else None), fmt(c[1] if c else None),
                    fmt(t[2] if t else None), fmt(c[2] if c else None), fmt(react),
                    fmt(t[3] if t else None), fmt(c[3] if c else None),
                    fmt(t[4] if t else None), fmt(c[4] if c else None)]) + "\n")
                n_rows += 1
                if t and c and t[1] >= args.min_cov and c[1] >= args.min_cov:
                    per_contig[key[0]].append((key[1], nt, t[1], c[1], t[2], c[2], react))

            for contig in contigs:
                rows = per_contig[contig]
                if not rows:
                    sum_out.write("\t".join(meta + [contig, "0"] + [""] * 13) + "\n")
                    continue
                react = [r[6] for r in rows]
                background = statistics.fmean(r[5] for r in rows)
                best = max(rows, key=lambda r: r[6])
                at_t = [r[6] for r in rows if r[1] == "T"]
                not_t = [r[6] for r in rows if r[1] != "T"]
                sum_out.write("\t".join(meta + [
                    contig, str(len(rows)),
                    fmt(statistics.median(r[2] for r in rows)),
                    fmt(statistics.median(r[3] for r in rows)),
                    fmt(background),
                    fmt(statistics.fmean(r[4] for r in rows)),
                    fmt(statistics.fmean(react)), fmt(best[6]), str(best[0]), best[1],
                    fmt(statistics.fmean(at_t) if at_t else None),
                    fmt(statistics.fmean(not_t) if not_t else None),
                    str(len(at_t)),
                    fmt(statistics.stdev(react) if len(react) > 1 else None),
                    str(sum(1 for x in react if x > 3 * background))]) + "\n")

    print(f"Manifest rows:      {len(manifest)}")
    print(f"Position rows:      {n_rows}")
    if n_empty:
        print(f"Manifest rows with no data in either sample: {n_empty} "
              f"(a sample with no reads on that reference)")
    print(f"Positions: {args.positions}")
    print(f"Summary:   {args.summary}")
    if n_rows == 0:
        sys.exit("Error: every file in the manifest was empty; there is nothing to "
                 "compare. Check the sweep's logs under logs/sweep/.")


if __name__ == "__main__":
    main()
