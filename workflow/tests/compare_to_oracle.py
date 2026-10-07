#!/usr/bin/env python3
"""
compare_to_oracle.py - Compare a perbase_signal output against the convert-TSV oracle.

The two paths are NOT expected to be byte-identical, and the reason is not a bug in
either: `uncalled4 convert` writes its TSV with `float_format="%.6g"`
(uncalled4 src/uncalled4/io/tsv.py:67), so the oracle's statistics were computed from
dtw.model_diff values already rounded to **6 significant digits**. The BAM-direct path
never serialises, so it uses the full float64 value. Where they differ, the BAM path
is the more accurate one.

Rounding to 6 significant digits is a RELATIVE error of <=5e-6, so it costs up to about
5e-6 absolute on a value near 1. Quantiles inherit the full error of one observation,
since each interpolates between two order statistics.

**Column 10 (mean squared deviation) needs a larger tolerance than the other five.**
d(x^2) = 2x*dx, so squaring amplifies each observation's rounding by 2|x| -- measured at
2x to 20x for |x| from 0.5 to 2.5. It therefore gets its own `--max-abs-sq`.

Coverage is the other half of it: averaging cancels these errors, so deep data agrees far
more closely than shallow. On js4022's full mtDNA sample (250-666x) column 10 differed by
at most 1e-6; on the nine edge-case BAMs (1-6 reads) it reaches 2.5e-5. The cases are a
deliberate worst case, so do not read a larger spread there as a regression.

For a verdict that needs no tolerance at all, use `--prove-precision`.

What MUST match exactly: row count, chromosome, position, nucleotide, coverage. Those
come from the decode, not the arithmetic, and any disagreement there is a real bug.

Usage:
    compare_to_oracle.py -a new.txt.gz -b oracle.txt.gz [--max-abs 1e-5]

    # prove the mechanism: re-run the decode, quantising each observation to 6
    # significant digits first, and expect byte-identity with the oracle
    compare_to_oracle.py --prove-precision -i uncalled4.bam -g ref.fa -s for \\
        -b oracle.txt.gz
"""

import argparse
import gzip
import os
import sys

import numpy as np

COLS = ["chr", "pos", "nt", "cov", "mean", "q25", "q75", "q025", "q975", "mean_sq"]
NUMERIC = list(range(4, 10))          # columns 5-10, 0-indexed


def load(path):
    opener = gzip.open if path.endswith(".gz") else open
    keys, vals = [], []
    with opener(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) != 10:
                sys.exit(f"Error: {path} has a row with {len(f)} fields, expected 10")
            keys.append((f[0], int(f[1]), f[2], int(f[3])))
            vals.append([float(x) for x in f[4:]])
    return keys, np.asarray(vals, dtype=np.float64)


def compare(path_a, path_b, max_abs, max_abs_sq):
    ka, va = load(path_a)
    kb, vb = load(path_b)

    print(f"A (new):    {path_a}   {len(ka)} rows")
    print(f"B (oracle): {path_b}   {len(kb)} rows")
    print()

    if len(ka) != len(kb):
        print(f"FAIL: row counts differ ({len(ka)} vs {len(kb)})")
        sa, sb = {k[:2] for k in ka}, {k[:2] for k in kb}
        only_a, only_b = sorted(sa - sb)[:5], sorted(sb - sa)[:5]
        if only_a:
            print(f"  only in A: {only_a}")
        if only_b:
            print(f"  only in B: {only_b}")
        return 1

    # --- the parts that must be exact ---------------------------------------
    bad = [i for i, (x, y) in enumerate(zip(ka, kb)) if x != y]
    if bad:
        print(f"FAIL: {len(bad)} rows differ in chr/pos/nt/coverage — this is a decode "
              f"bug, not rounding.")
        for i in bad[:5]:
            print(f"  row {i}: A={ka[i]}  B={kb[i]}")
        return 1
    print("chr / pos / nt / coverage: IDENTICAL on all rows "
          "(decode, strand filter, coordinates and NaN handling all agree)")
    print()

    # --- the parts that may differ by rounding -------------------------------
    print(f"{'column':<10} {'exact':>8} {'differ':>8} {'max|diff|':>12} "
          f"{'at pos':>9} {'tol':>9} {'':>4}")
    failed = []
    for j, name in zip(NUMERIC, COLS[4:]):
        d = np.abs(va[:, j - 4] - vb[:, j - 4])
        n_diff = int((d > 0).sum())
        k = int(d.argmax())
        # Column 10 is mean(x^2); squaring amplifies each observation's rounding
        # by 2|x|, so it gets its own bound. See the module docstring.
        tol = max_abs_sq if name == "mean_sq" else max_abs
        ok = float(d[k]) <= tol
        if not ok:
            failed.append(name)
        print(f"{name:<10} {len(d) - n_diff:>8} {n_diff:>8} {d[k]:>12.3e} "
              f"{ka[k][1]:>9} {tol:>9.1e} {'' if ok else '  OVER'}")

    print()
    if not failed:
        print(f"PASS: every column within tolerance (linear {max_abs:.1e}, "
              f"squared {max_abs_sq:.1e}) — consistent with the oracle's "
              f"6-significant-digit TSV quantisation.")
        return 0
    print(f"FAIL: over tolerance in {', '.join(failed)}. Larger than %.6g rounding "
          f"explains. Run --prove-precision before widening any tolerance: if that "
          f"reports byte-identity, the cause is rounding after all and the tolerance "
          f"is what needs revisiting.")
    return 1


def prove_precision(bam, genome, strand, oracle, script):
    """Re-run the decode with each observation rounded to 6 significant digits.

    If that reproduces the oracle byte-for-byte, the remaining differences are
    entirely the TSV's float_format and nothing else.
    """
    import importlib.util
    spec = importlib.util.spec_from_file_location("psd", script)
    psd = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(psd)

    import pysam
    fasta = pysam.FastaFile(genome)
    seq_cache = {}

    out_path = oracle + ".requantised.txt.gz"
    written = 0
    with gzip.open(out_path, "wt") as out:
        def emit(ref, _w, pos, val):
            nonlocal written
            if ref not in seq_cache:
                seq_cache.clear()
                seq_cache[ref] = fasta.fetch(ref).upper()
            for line in psd.format_window(ref, pos, val, seq_cache[ref], 1):
                out.write(line)
                written += 1

        acc = psd.WindowAccumulator(10000, emit)
        for ref, start, pos, val in psd.iter_read_observations(
                bam, genome, strand == "rev"):
            # exactly what uncalled4's TSV writer does to every value
            val = np.array([float("%.6g" % v) for v in val], dtype=np.float64)
            acc.set_ref(ref)
            acc.add(pos, val, start)
        acc.close()
    fasta.close()

    a = gzip.open(out_path, "rt").read()
    b = gzip.open(oracle, "rt").read()
    print(f"re-quantised output: {out_path} ({written} rows)")
    if a == b:
        print("PROVEN: byte-identical to the oracle once each observation is rounded "
              "to 6 significant digits. The only difference between the two paths is "
              "the TSV's float_format, and the BAM path is the more accurate one.")
        return 0
    la, lb = a.splitlines(), b.splitlines()
    print(f"NOT identical: {sum(1 for x, y in zip(la, lb) if x != y)} of {len(la)} "
          f"rows still differ — something beyond TSV rounding is going on.")
    for i, (x, y) in enumerate(zip(la, lb)):
        if x != y:
            print(f"  first diff row {i}:\n    requantised {x}\n    oracle      {y}")
            break
    return 1


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-a", "--new", help="Output from the BAM-direct path")
    p.add_argument("-b", "--oracle", required=True, help="Output from the convert-TSV path")
    p.add_argument("--max-abs", type=float, default=1e-5,
                   help="Largest absolute difference to accept in the mean and the four "
                        "quantiles (default: 1e-5)")
    p.add_argument("--max-abs-sq", type=float, default=5e-5,
                   help="Largest absolute difference to accept in mean_sq (default: "
                        "5e-5). Higher than --max-abs because d(x^2)=2x*dx amplifies "
                        "the oracle's %%.6g rounding by 2|x|; shallow coverage removes "
                        "the averaging that would otherwise cancel it.")
    p.add_argument("--prove-precision", action="store_true",
                   help="Re-run the decode with 6-significant-digit rounding and "
                        "require byte-identity with the oracle")
    p.add_argument("-i", "--bam", help="With --prove-precision: the Uncalled4 BAM")
    p.add_argument("-g", "--genome", help="With --prove-precision: reference FASTA")
    p.add_argument("-s", "--strand", choices=["for", "rev"],
                   help="With --prove-precision: strand")
    p.add_argument("--script",
                   default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                        "..", "scripts",
                                        "perbase_signal_deviation.py"),
                   help="Path to perbase_signal_deviation.py")
    args = p.parse_args()

    if args.prove_precision:
        missing = [f for f, v in (("-i", args.bam), ("-g", args.genome),
                                  ("-s", args.strand)) if not v]
        if missing:
            p.error(f"--prove-precision needs {', '.join(missing)}")
        sys.exit(prove_precision(args.bam, args.genome, args.strand,
                                 args.oracle, args.script))

    if not args.new:
        p.error("-a/--new is required without --prove-precision")
    sys.exit(compare(args.new, args.oracle, args.max_abs, args.max_abs_sq))


if __name__ == "__main__":
    main()
