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
  Which positions get a row is set by --position-rows and --regions (below).

summary (TSV) - one row per manifest row x contig. Over positions where both samples
have at least --min-cov reads:
    n_positions
    cov_treated_median, cov_control_median
    background            mean frac_control - the false-positive rate actually
                          realised on the held-out control
    frac_treated_mean
    reactivity_mean, reactivity_max, pos_of_max, nt_of_max
    reactivity_mean_T, reactivity_mean_nonT, n_T
                          split by the base on the read's strand: T on 'for', and
                          on 'rev' an A in the nt column (nt is always the forward
                          reference base, so a T on the reverse strand reads A).
                          Permanganate reacts mainly with T, but one modified T
                          shifts the current at its neighbours too, so non-T is
                          not a clean negative.
    reactivity_sd         spread of reactivity across positions
    n_reactive            positions with reactivity > 3 x background
  Then, pooled over every position both samples scored, whatever its coverage
  (observation-weighted, so usable at depths where few positions reach --min-cov):
    n_shared              positions scored in both samples
    obs_treated, obs_control
                          observations (read x position) in those positions
    pooled_background     modified observations / observations, control
    pooled_frac_treated   the same, treated
    pooled_reactivity     pooled_frac_treated - pooled_background
    pooled_reactivity_T, pooled_reactivity_nonT
    pooled_mean_z_treated, pooled_mean_z_control
    pooled_mean_absz_treated, pooled_mean_absz_control
                          threshold-free: the shift and the size of the signal change
    pooled_absz_excess_T, pooled_absz_excess_nonT
                          treated - control mean |z|, at T and at non-T
    scope                 'all', or 'regions' for the extra row per contig written
                          over the --regions intervals only
  The pooled columns ignore that observations are not independent (neighbouring
  k-mers share a base; reads cover many positions), so no p-value is given.

Nothing here chooses between parameter sets; it lays them out to be compared.

Rows of the positions table:
  --position-rows all     every position (what model oligos want)
  --position-rows none    header only
  --position-rows auto    (default) all, unless a treated file has more than
                          AUTO_MAX_POSITIONS positions - a genome would give one row
                          per position per manifest row, billions of rows - then none
  --regions BED           only positions inside these intervals (BED: 0-based start,
                          end exclusive; extra columns ignored), and adds the
                          scope=regions summary rows. Overrides auto and all;
                          with --position-rows none the regions get their summary
                          rows but no position rows (sites too thinly covered to
                          read one by one).

Each manifest row is independent, so they are spread over --threads processes. A
call file is held as numpy arrays (~40 bytes per position), and each process keeps
only the last two it read.

Usage:
    mod_sweep_table.py -m manifest.tsv -p positions.tsv.gz -s summary.tsv \
        [-c 10] [-t 8] [--position-rows auto|all|none] [--regions sites.bed]
"""

import os
import sys
import csv
import gzip
import math
import argparse
from array import array
from collections import OrderedDict
from multiprocessing import Pool

import numpy as np

AUTO_MAX_POSITIONS = 200000

POS_COLS = ["chr", "pos", "nt", "cov_treated", "cov_control", "frac_treated",
            "frac_control", "reactivity", "mean_z_treated", "mean_z_control",
            "mean_absz_treated", "mean_absz_control"]
SUM_COLS = ["chr", "n_positions", "cov_treated_median", "cov_control_median",
            "background", "frac_treated_mean", "reactivity_mean", "reactivity_max",
            "pos_of_max", "nt_of_max", "reactivity_mean_T", "reactivity_mean_nonT",
            "n_T", "reactivity_sd", "n_reactive"]
POOLED_COLS = ["n_shared", "obs_treated", "obs_control", "pooled_background",
               "pooled_frac_treated", "pooled_reactivity", "pooled_reactivity_T",
               "pooled_reactivity_nonT", "pooled_mean_z_treated",
               "pooled_mean_z_control", "pooled_mean_absz_treated",
               "pooled_mean_absz_control", "pooled_absz_excess_T",
               "pooled_absz_excess_nonT", "scope"]


def parse_args():
    p = argparse.ArgumentParser(
        description="Tabulate a modification-calling parameter sweep",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    p.add_argument("-m", "--manifest", required=True)
    p.add_argument("-p", "--positions", required=True, help="Output (.tsv.gz)")
    p.add_argument("-s", "--summary", required=True, help="Output (.tsv)")
    p.add_argument("-c", "--min-cov", type=int, default=10,
                   help="Coverage both samples need for a position to enter the "
                        "per-position summary columns (default: 10)")
    p.add_argument("-t", "--threads", type=int, default=1,
                   help="Processes (default: 1)")
    p.add_argument("--position-rows", choices=["auto", "all", "none"], default="auto",
                   help="Which positions get a row in the positions table "
                        "(default: auto)")
    p.add_argument("--regions", default=None,
                   help="BED of intervals: positions table and scope=regions "
                        "summary rows restricted to these")
    return p.parse_args()


class Calls:
    """One perbase_mod_calls.py output, as numpy arrays per contig sorted by position."""

    def __init__(self, path):
        self.order = []
        self.contigs = {}
        # Typed columns, not a list of tuples: a 4.2 Mb strand is ~200 MB this way
        # and ~1.7 GB as Python objects.
        cols = {}
        with gzip.open(path, "rt") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                c = cols.get(f[0])
                if c is None:
                    self.order.append(f[0])
                    c = cols[f[0]] = (bytearray(), array("q"), array("q"), array("d"),
                                      array("q"), array("d"), array("d"))
                c[0].extend(f[2].encode())
                c[1].append(int(f[1]))
                c[2].append(int(f[3]))
                c[3].append(float(f[4]))
                c[4].append(int(f[5]))
                c[5].append(float(f[6]))
                c[6].append(float(f[7]))
        for contig, (nt, pos, cov, frac, n_mod, z, absz) in cols.items():
            arrays = dict(
                nt=np.frombuffer(bytes(nt), dtype="S1"),
                pos=np.frombuffer(pos, dtype=np.int64), cov=np.frombuffer(cov, dtype=np.int64),
                frac=np.frombuffer(frac), n_mod=np.frombuffer(n_mod, dtype=np.int64),
                z=np.frombuffer(z), absz=np.frombuffer(absz))
            step = np.diff(arrays["pos"])
            if np.any(step < 0):
                # perbase_mod_calls.py writes in position order; sort only if not.
                idx = np.argsort(arrays["pos"], kind="stable")
                arrays = {k: v[idx] for k, v in arrays.items()}
                step = np.diff(arrays["pos"])
            if np.any(step == 0):
                sys.exit(f"Error: {path} has a position twice on {contig}.")
            self.contigs[contig] = arrays
        self.n_positions = sum(a["pos"].size for a in self.contigs.values())


def fmt(x):
    return "" if x is None else (f"{x:.6f}" if isinstance(x, float) else str(x))


def py_median(values):
    """statistics.median on ints, which this table has always used: the middle int
    for an odd count, the float mean of the middle two for an even one."""
    v = np.sort(values)
    n = v.size
    if n % 2:
        return int(v[n // 2])
    return (int(v[n // 2 - 1]) + int(v[n // 2])) / 2


def fmean(values):
    return math.fsum(values.tolist()) / values.size


def read_regions(path):
    regions = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            if len(f) < 3:
                sys.exit(f"Error: {path}: a BED line needs 3 tab-separated fields: {line!r}")
            regions.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    merged = {}
    for contig, iv in regions.items():
        iv.sort()
        out = [list(iv[0])]
        for s, e in iv[1:]:
            if s <= out[-1][1]:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e])
        merged[contig] = (np.array([s for s, _ in out], dtype=np.int64),
                          np.array([e for _, e in out], dtype=np.int64))
    return merged


def in_regions(contig, pos, regions):
    """1-based positions -> mask of those inside the (0-based, end-exclusive) BED."""
    if contig not in regions:
        return np.zeros(pos.size, dtype=bool)
    starts, ends = regions[contig]
    i = np.searchsorted(starts, pos - 1, side="right") - 1
    return (i >= 0) & (pos - 1 < ends[np.maximum(i, 0)])


# Per-process state: options, regions, and a two-file cache.
STATE = {}


def init_worker(min_cov, rows_mode, regions):
    STATE.update(min_cov=min_cov, rows_mode=rows_mode, regions=regions,
                 cache=OrderedDict())


def load(path):
    cache = STATE["cache"]
    if path in cache:
        cache.move_to_end(path)
        return cache[path]
    if not os.path.exists(path):
        raise FileNotFoundError(f"{path} is listed in the manifest but does not exist.")
    calls = Calls(path)
    cache[path] = calls
    while len(cache) > 2:
        cache.popitem(last=False)
    return calls


def empty_side(n):
    return dict(nt=np.full(n, b"N", dtype="S1"), pos=np.zeros(n, dtype=np.int64),
                cov=np.zeros(n, dtype=np.int64), frac=np.zeros(n),
                n_mod=np.zeros(n, dtype=np.int64), z=np.zeros(n), absz=np.zeros(n))


def join(t, c):
    """Union of positions of one contig: (pos, has_t, has_c, t arrays, c arrays),
    each side's arrays aligned to the union (zeros where absent)."""
    t_pos = t["pos"] if t else np.zeros(0, dtype=np.int64)
    c_pos = c["pos"] if c else np.zeros(0, dtype=np.int64)
    pos = np.union1d(t_pos, c_pos)
    sides = []
    for side, side_pos in ((t, t_pos), (c, c_pos)):
        has = np.isin(pos, side_pos, assume_unique=True)
        aligned = empty_side(pos.size)
        if side:
            for k, v in side.items():
                aligned[k][has] = v
        sides.append((has, aligned))
    return pos, sides[0][0], sides[1][0], sides[0][1], sides[1][1]


def pooled(mask, is_t, ta, ca):
    """The observation-weighted columns over the positions in mask."""
    n = int(mask.sum())
    if n == 0:
        return [str(n)] + [""] * 13

    def frac(side, m):
        obs = int(side["cov"][m].sum())
        return int(side["n_mod"][m].sum()) / obs if obs else None

    def weighted(side, key, m):
        obs = int(side["cov"][m].sum())
        return math.fsum((side[key][m] * side["cov"][m]).tolist()) / obs if obs else None

    def diff(a, b):
        return None if a is None or b is None else a - b

    mt, mn = mask & is_t, mask & ~is_t
    bg, ft = frac(ca, mask), frac(ta, mask)
    return [str(n), str(int(ta["cov"][mask].sum())), str(int(ca["cov"][mask].sum())),
            fmt(bg), fmt(ft), fmt(diff(ft, bg)),
            fmt(diff(frac(ta, mt), frac(ca, mt))), fmt(diff(frac(ta, mn), frac(ca, mn))),
            fmt(weighted(ta, "z", mask)), fmt(weighted(ca, "z", mask)),
            fmt(weighted(ta, "absz", mask)), fmt(weighted(ca, "absz", mask)),
            fmt(diff(weighted(ta, "absz", mt), weighted(ca, "absz", mt))),
            fmt(diff(weighted(ta, "absz", mn), weighted(ca, "absz", mn)))]


def per_position(mask, pos, nt, is_t, ta, ca):
    """The per-position (min-cov) columns over the positions in mask."""
    n = int(mask.sum())
    if n == 0:
        return [str(n)] + [""] * 13
    react = ta["frac"][mask] - ca["frac"][mask]
    background = fmean(ca["frac"][mask])
    best = int(np.argmax(react))   # first of equal maxima, as max() picked
    at_t, not_t = react[is_t[mask]], react[~is_t[mask]]
    if react.size > 1:
        mean = fmean(react)
        sd = math.sqrt(math.fsum(((react - mean) ** 2).tolist()) / (react.size - 1))
    else:
        sd = None
    return [str(n), fmt(py_median(ta["cov"][mask])), fmt(py_median(ca["cov"][mask])),
            fmt(background), fmt(fmean(ta["frac"][mask])), fmt(fmean(react)),
            fmt(float(react[best])), str(int(pos[mask][best])),
            nt[mask][best].decode(),
            fmt(fmean(at_t) if at_t.size else None),
            fmt(fmean(not_t) if not_t.size else None), str(int(at_t.size)), fmt(sd),
            str(int(np.sum(react > 3 * background)))]


def process_unit(args):
    unit, meta = args
    min_cov, rows_mode, regions = STATE["min_cov"], STATE["rows_mode"], STATE["regions"]
    treated = load(unit["treatment_file"])
    control = load(unit["control_file"])
    contigs = treated.order + [c for c in control.order if c not in treated.order]
    pos_lines, sum_lines = [], []
    n_rows = 0
    for contig in contigs:
        pos, has_t, has_c, ta, ca = join(treated.contigs.get(contig),
                                          control.contigs.get(contig))
        nt = np.where(has_t, ta["nt"], ca["nt"])
        # nt is the forward reference base on both strands; a T on the reverse
        # strand is an A there (js4037: testing for T on 'rev' swapped the split).
        is_t = nt == (b"A" if unit.get("strand") == "rev" else b"T")
        shared = has_t & has_c
        enough = shared & (ta["cov"] >= min_cov) & (ca["cov"] >= min_cov)

        if rows_mode == "none":
            write = np.zeros(pos.size, dtype=bool)
        elif regions is not None:
            write = in_regions(contig, pos, regions)
        else:
            write = np.full(pos.size, rows_mode == "all")
        for k in np.flatnonzero(write):
            t, c = bool(has_t[k]), bool(has_c[k])
            pos_lines.append("\t".join(meta + [
                contig, str(int(pos[k])), nt[k].decode(),
                str(int(ta["cov"][k])) if t else "", str(int(ca["cov"][k])) if c else "",
                fmt(float(ta["frac"][k])) if t else "",
                fmt(float(ca["frac"][k])) if c else "",
                fmt(float(ta["frac"][k] - ca["frac"][k])) if t and c else "",
                fmt(float(ta["z"][k])) if t else "", fmt(float(ca["z"][k])) if c else "",
                fmt(float(ta["absz"][k])) if t else "",
                fmt(float(ca["absz"][k])) if c else ""]))
        n_rows += int(write.sum())

        scopes = [("all", np.ones(pos.size, dtype=bool))]
        if regions is not None:
            scopes.append(("regions", in_regions(contig, pos, regions)))
        for scope, sel in scopes:
            sum_lines.append("\t".join(
                meta + [contig] + per_position(enough & sel, pos, nt, is_t, ta, ca)
                + pooled(shared & sel, is_t, ta, ca) + [scope]))
    empty = not treated.contigs and not control.contigs
    return pos_lines, sum_lines, n_rows, empty


def main():
    args = parse_args()
    with open(args.manifest) as fh:
        manifest = list(csv.DictReader(fh, delimiter="\t"))
    if not manifest:
        sys.exit(f"Error: manifest {args.manifest} has no rows.")
    meta_cols = [c for c in manifest[0] if c not in ("treatment_file", "control_file")]
    for unit in manifest:
        for key in ("treatment_file", "control_file"):
            if not os.path.exists(unit[key]):
                sys.exit(f"Error: {unit[key]} is listed in the manifest but does not exist.")

    regions = read_regions(args.regions) if args.regions else None
    rows_mode = args.position_rows
    if rows_mode == "none":
        print("Position rows:      none")
    elif regions is not None:
        print(f"Position rows:      inside {sum(s.size for s, _ in regions.values())} "
              f"intervals of {args.regions}")
    elif rows_mode == "auto":
        # One treated file per genome decides: the references are what sets the size.
        biggest = 0
        for genome in dict.fromkeys(u.get("genome", "") for u in manifest):
            first = next(u for u in manifest if u.get("genome", "") == genome)
            biggest = max(biggest, Calls(first["treatment_file"]).n_positions)
        rows_mode = "all" if biggest <= AUTO_MAX_POSITIONS else "none"
        print(f"Position rows:      {rows_mode} (auto: largest treated file has "
              f"{biggest} positions; limit {AUTO_MAX_POSITIONS})")
    else:
        print(f"Position rows:      {rows_mode}")

    jobs = [(unit, [unit[c] for c in meta_cols]) for unit in manifest]
    n_rows = n_empty = 0
    with gzip.open(args.positions, "wt") as pos_out, open(args.summary, "w") as sum_out:
        pos_out.write("\t".join(meta_cols + POS_COLS) + "\n")
        sum_out.write("\t".join(meta_cols + SUM_COLS + POOLED_COLS) + "\n")
        initargs = (args.min_cov, rows_mode, regions)
        if args.threads > 1:
            # Contiguous chunks keep neighbouring manifest rows - which share their
            # control file - in one process's cache.
            chunk = max(1, len(jobs) // (args.threads * 4))
            pool = Pool(args.threads, initializer=init_worker, initargs=initargs)
            results = pool.imap(process_unit, jobs, chunksize=chunk)
        else:
            pool = None
            init_worker(*initargs)
            results = map(process_unit, jobs)
        try:
            for pos_lines, sum_lines, rows, empty in results:
                for line in pos_lines:
                    pos_out.write(line + "\n")
                for line in sum_lines:
                    sum_out.write(line + "\n")
                n_rows += rows
                n_empty += empty
        finally:
            if pool is not None:
                pool.close()
                pool.join()

    print(f"Manifest rows:      {len(manifest)}")
    print(f"Position rows:      {n_rows}")
    if n_empty:
        print(f"Manifest rows with no data in either sample: {n_empty} "
              f"(a sample with no reads on that reference)")
    print(f"Positions: {args.positions}")
    print(f"Summary:   {args.summary}")
    if n_empty == len(manifest):
        sys.exit("Error: every file in the manifest was empty; there is nothing to "
                 "compare. Check the sweep's logs under logs/sweep/.")


if __name__ == "__main__":
    main()
