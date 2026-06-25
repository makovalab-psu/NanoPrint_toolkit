#!/usr/bin/env python3
"""
perbase_error.py - Compute per-base error statistics from samtools mpileup output.

Reads mpileup from stdin; writes 9-column tab-delimited output to stdout.

Input: samtools mpileup -f <ref> -B -Q 0 <bam>
    mpileup columns used: 1=CHROM, 2=POS (1-based), 3=REF, 6=base qualities

Output columns (tab-delimited):
    1. Chromosome
    2. Position (1-based)
    3. Reference nucleotide
    4. Coverage
    5. Mean error probability  (mean of per-read 10^(-Q/10))
    6. Q25  — 0.25 quantile (lower 50% CI bound)
    7. Q75  — 0.75 quantile (upper 50% CI bound)
    8. Q025 — 0.025 quantile (lower 95% CI bound)
    9. Q975 — 0.975 quantile (upper 95% CI bound)
"""

import sys


def _quantile(sorted_vals, q):
    n = len(sorted_vals)
    idx = q * (n - 1)
    lo = int(idx)
    hi = lo + 1
    if hi >= n:
        return sorted_vals[-1]
    return sorted_vals[lo] + (idx - lo) * (sorted_vals[hi] - sorted_vals[lo])


def main():
    for line in sys.stdin:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 6:
            continue
        chrom, pos, ref, quals = fields[0], fields[1], fields[2], fields[5]
        if not quals:
            continue

        probs = sorted(10 ** (-(ord(ch) - 33) / 10) for ch in quals)
        n = len(probs)

        mean_err = sum(probs) / n
        q25  = _quantile(probs, 0.25)
        q75  = _quantile(probs, 0.75)
        q025 = _quantile(probs, 0.025)
        q975 = _quantile(probs, 0.975)

        print(
            f"{chrom}\t{pos}\t{ref}\t{n}\t{mean_err:.6f}"
            f"\t{q25:.6f}\t{q75:.6f}\t{q025:.6f}\t{q975:.6f}"
        )


if __name__ == "__main__":
    main()
