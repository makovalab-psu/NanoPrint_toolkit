#!/usr/bin/env python3
"""
read_summary_table.py - Are modifications concentrated in some reads?

Compares the per-read call rates of a treated sample with those of its held-out
control, from perbase_mod_calls.py --read-summary files. A per-position table
cannot tell "every molecule modified at 3%" from "10% of molecules at 30%"; this
can. js4037: permanganate is expected to reach some cells (damaged envelopes) and
not others, so treated reads should carry a tail of high call rates that the
control does not.

Input: a manifest (TSV with a header) with at least the columns
    set  relationship  strand  treatment_reads  control_reads
and any others, which are copied to the output. Written by mod_sweep.smk
(rule sweep_read_table).

A read's call rate is n_called / n_scored. Only reads with at least --min-scored
scored positions are used: a short read's rate is mostly counting noise.

The control's false calls are not binomial - low-quality reads have more - so the
treated reads are judged against the control's own distribution, within bins of
the read's mean basecall quality (qs). Bin edges are the control's qs quantiles
(--qs-bins bins; one bin when reads have no qs). In each bin, a read is
'reactive' when its call rate is above the control's --tail quantile there. The
control is then reactive at 1 - tail by construction, and everything treated above
that is the excess.

summary (TSV) - one row per manifest row:
    the manifest's columns except the two file paths, then for treated (_t) and
    control (_c):
      n_reads             reads used (>= --min-scored positions)
      rate_mean, rate_median, rate_p90, rate_p99
                          per-read call rate
      dispersion          variance of the per-read rate / the variance a binomial
                          would give each read at the sample's pooled rate.
                          1 = calls spread evenly over reads; much more than the
                          control's = some reads carry more than their share
    frac_reactive_t       treated reads above the control's tail cutoff
    excess_reactive       frac_reactive_t - (1 - tail)
    calls_in_reactive_t   share of all treated calls that sit in reactive reads
    rate_T_reactive, rate_nonT_reactive
    rate_T_rest, rate_nonT_rest
    rate_T_c, rate_nonT_c
                          pooled call rate at T and non-T positions, in the
                          treated reactive reads, the other treated reads, and the
                          control. Permanganate oxidises T: a real tail should be
                          T-biased; a tail with T ~ non-T is damage or junk.

hist (TSV, long) - the per-read rate distributions, for plotting:
    the manifest's columns except the file paths, sample (treated|control),
    rate_bin (lower edge), n_reads
"""

import sys
import csv
import gzip
import argparse

import numpy as np

FILE_COLUMNS = ("treatment_reads", "control_reads")


def parse_args():
    p = argparse.ArgumentParser(
        description="Compare per-read modification call rates, treated vs control")
    p.add_argument("-m", "--manifest", required=True)
    p.add_argument("-s", "--summary", required=True, help="Summary TSV to write")
    p.add_argument("--hist", required=True, help="Histogram TSV to write")
    p.add_argument("--min-scored", type=int, default=1000,
                   help="Reads with fewer scored positions are not used (default: 1000)")
    p.add_argument("--qs-bins", type=int, default=5,
                   help="Quality bins the control's cutoff is taken within (default: 5)")
    p.add_argument("--tail", type=float, default=0.99,
                   help="Control quantile a reactive read must exceed (default: 0.99)")
    p.add_argument("--bin-width", type=float, default=0.005,
                   help="Width of the call-rate histogram bins (default: 0.005)")
    args = p.parse_args()
    if not 0 < args.tail < 1:
        p.error("--tail must be between 0 and 1")
    return args


def load_reads(path, min_scored):
    """Columns of a --read-summary file as arrays, reads below min_scored dropped."""
    cols = {k: [] for k in ("qs", "n_scored", "n_called", "n_T", "n_T_called")}
    with gzip.open(path, "rt") if path.endswith(".gz") else open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cols["qs"].append(np.nan if row["qs"] == "NA" else float(row["qs"]))
            for k in ("n_scored", "n_called", "n_T", "n_T_called"):
                cols[k].append(int(row[k]))
    out = {k: np.asarray(v, dtype=np.float64) for k, v in cols.items()}
    keep = out["n_scored"] >= min_scored
    out = {k: v[keep] for k, v in out.items()}
    out["rate"] = (out["n_called"] / out["n_scored"]) if keep.any() else np.empty(0)
    return out


def qs_bin_of(qs, edges):
    """Bin index of each read; reads without qs go to the middle bin."""
    if edges is None:
        return np.zeros(qs.size, dtype=np.int64)
    b = np.searchsorted(edges, qs, side="right") - 1
    b = np.clip(b, 0, edges.size - 2)
    b[np.isnan(qs)] = (edges.size - 2) // 2
    return b


def rate_stats(r):
    if r["rate"].size == 0:
        return dict(n_reads=0, rate_mean="NA", rate_median="NA", rate_p90="NA",
                    rate_p99="NA", dispersion="NA")
    rate, n = r["rate"], r["n_scored"]
    pooled = r["n_called"].sum() / n.sum()
    binom_var = np.mean(pooled * (1 - pooled) / n)
    return dict(n_reads=rate.size,
                rate_mean=f"{rate.mean():.5f}",
                rate_median=f"{np.median(rate):.5f}",
                rate_p90=f"{np.quantile(rate, 0.9):.5f}",
                rate_p99=f"{np.quantile(rate, 0.99):.5f}",
                dispersion=f"{rate.var() / binom_var:.3f}" if binom_var > 0 else "NA")


def t_rates(r, mask):
    """Pooled call rate at T and at non-T positions over the reads in mask."""
    n_t, c_t = r["n_T"][mask].sum(), r["n_T_called"][mask].sum()
    n_o = r["n_scored"][mask].sum() - n_t
    c_o = r["n_called"][mask].sum() - c_t
    fmt = lambda c, n: f"{c / n:.5f}" if n > 0 else "NA"
    return fmt(c_t, n_t), fmt(c_o, n_o)


def main():
    args = parse_args()
    with open(args.manifest) as fh:
        units = list(csv.DictReader(fh, delimiter="\t"))
    if not units:
        sys.exit(f"Error: manifest {args.manifest} has no rows.")
    missing = [c for c in ("set", "relationship", "strand") + FILE_COLUMNS
               if c not in units[0]]
    if missing:
        sys.exit(f"Error: manifest {args.manifest} lacks column(s): {', '.join(missing)}")
    keys = [c for c in units[0] if c not in FILE_COLUMNS]

    stat_cols = ["n_reads", "rate_mean", "rate_median", "rate_p90", "rate_p99",
                 "dispersion"]
    summary_cols = (keys + [f"{c}_t" for c in stat_cols] + [f"{c}_c" for c in stat_cols]
                    + ["frac_reactive_t", "excess_reactive", "calls_in_reactive_t",
                       "rate_T_reactive", "rate_nonT_reactive", "rate_T_rest",
                       "rate_nonT_rest", "rate_T_c", "rate_nonT_c"])
    edges_hist = np.arange(0, 1 + args.bin_width, args.bin_width)

    print("=== Per-read Call Rates ===")
    print(f"Manifest:      {args.manifest} ({len(units)} rows)")
    print(f"Reads used:    >= {args.min_scored} scored positions")
    print(f"Reactive:      above the control's {args.tail} quantile, within "
          f"{args.qs_bins} qs bins")
    print()

    with open(args.summary, "w") as out_s, open(args.hist, "w") as out_h:
        out_s.write("\t".join(summary_cols) + "\n")
        out_h.write("\t".join(keys + ["sample", "rate_bin", "n_reads"]) + "\n")
        for unit in units:
            t = load_reads(unit["treatment_reads"], args.min_scored)
            c = load_reads(unit["control_reads"], args.min_scored)
            row = {k: unit[k] for k in keys}
            row.update({f"{k}_t": v for k, v in rate_stats(t).items()})
            row.update({f"{k}_c": v for k, v in rate_stats(c).items()})

            if t["rate"].size and c["rate"].size:
                qs_ok = c["qs"][~np.isnan(c["qs"])]
                edges = None
                if args.qs_bins > 1 and qs_ok.size >= 10 * args.qs_bins:
                    edges = np.unique(np.quantile(qs_ok, np.linspace(0, 1, args.qs_bins + 1)))
                    edges = edges if edges.size >= 3 else None
                bins_c, bins_t = qs_bin_of(c["qs"], edges), qs_bin_of(t["qs"], edges)
                reactive = np.zeros(t["rate"].size, dtype=bool)
                for b in np.unique(bins_t):
                    in_c = c["rate"][bins_c == b]
                    if in_c.size == 0:
                        in_c = c["rate"]
                    reactive[bins_t == b] = t["rate"][bins_t == b] > np.quantile(in_c, args.tail)
                frac = reactive.mean()
                total = t["n_called"].sum()
                row.update(frac_reactive_t=f"{frac:.5f}",
                           excess_reactive=f"{frac - (1 - args.tail):.5f}",
                           calls_in_reactive_t=(f"{t['n_called'][reactive].sum() / total:.5f}"
                                                if total else "NA"))
                row["rate_T_reactive"], row["rate_nonT_reactive"] = t_rates(t, reactive)
                row["rate_T_rest"], row["rate_nonT_rest"] = t_rates(t, ~reactive)
            row["rate_T_c"], row["rate_nonT_c"] = t_rates(c, np.ones(c["rate"].size, bool))
            out_s.write("\t".join(str(row.get(k, "NA")) for k in summary_cols) + "\n")

            for sample, r in (("treated", t), ("control", c)):
                counts, _ = np.histogram(np.minimum(r["rate"], 1.0), bins=edges_hist)
                for lo, n in zip(edges_hist[:-1], counts):
                    if n:
                        out_h.write("\t".join([unit[k] for k in keys]
                                              + [sample, f"{lo:.4f}", str(n)]) + "\n")
            print(f"{unit['set']}  {unit['relationship']}  {unit['strand']}: "
                  f"{row['n_reads_t']} treated / {row['n_reads_c']} control reads, "
                  f"reactive {row.get('frac_reactive_t', 'NA')}")
    print("\nDone.")


if __name__ == "__main__":
    main()
