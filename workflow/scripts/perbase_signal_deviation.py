#!/usr/bin/env python3
"""
perbase_signal_deviation.py - Compute per-base pore model signal deviation
from an Uncalled4 DTW alignment TSV.

The TSV is produced by Uncalled4_convert_tsv.sh using:
  --tsv-cols "dtw.current,dtw.current_sd,dtw.start,dtw.length,dtw.model_diff"
Note: dtw.base is NOT requested — it is not a valid layer in current uncalled4 versions.

Per js4004/workflow/scripts/parse_uncalled4_output.py:
  - Column names may contain dots (e.g. "dtw.model_diff") — normalize to underscores
  - The chromosome column may be named "ref", "seq_name", or "chr"
  - The position column may be named "pos", "seq_pos", or "ref_pos"
  - dtw_model_diff = model current - observed current (pA); per official uncalled4 docs,
    this is predicted MINUS observed (positive = observed current lower than pore model)

Nucleotide source (in priority order):
  1. dtw.base column from uncalled4 TSV (if present and not all-NaN)
     - May be a letter (A/C/G/T) or an integer (0/1/2/3) — both are handled
  2. pysam.FastaFile lookup from reference FASTA (fallback)
     - Only loaded if dtw.base is absent or unusable

Output format (tab-delimited, gzipped):
    1. Chromosome
    2. Position (1-based)
    3. Nucleotide
    4. Coverage (reads contributing to this position)
    5. Mean signal deviation (mean dtw.model_diff, pA)
    6. Q25 — 0.25 quantile of dtw.model_diff (lower 50% CI bound, pA)
    7. Q75 — 0.75 quantile of dtw.model_diff (upper 50% CI bound, pA)
    8. Q025 — 0.025 quantile of dtw.model_diff (lower 95% CI bound, pA)
    9. Q975 — 0.975 quantile of dtw.model_diff (upper 95% CI bound, pA)

Usage:
    perbase_signal_deviation.py -i <dtw.tsv> -g <genome.fa> -o <output.txt.gz>
"""

import sys
import gzip
import argparse

try:
    import pandas as pd
except ImportError:
    sys.exit("Error: pandas not found. Install with: conda install pandas")


# Integer encoding used by uncalled4's "binarized" dtw.base column (if applicable)
_INT_TO_BASE = {0: "A", 1: "C", 2: "G", 3: "T"}


def parse_args():
    p = argparse.ArgumentParser(
        description="Per-base signal deviation from Uncalled4 DTW TSV",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    p.add_argument("-i", "--tsv", required=True,
                   help="Input Uncalled4 DTW TSV (strand-filtered, all chromosomes)")
    p.add_argument("-g", "--genome", required=True,
                   help="Reference genome FASTA (fallback for nucleotide lookup if dtw.base absent)")
    p.add_argument("-o", "--output", required=True,
                   help="Output file (.txt.gz)")
    p.add_argument("-c", "--min-cov", type=int, default=1,
                   help="Minimum coverage to emit a position (default: 1)")
    return p.parse_args()


def normalize_columns(df):
    """Rename dot-separated Uncalled4 columns and handle ref/pos name variants."""
    df = df.rename(columns={c: c.replace(".", "_") for c in df.columns})

    if "ref" not in df.columns:
        for alt in ("seq_name", "chr"):
            if alt in df.columns:
                df = df.rename(columns={alt: "ref"})
                break

    if "pos" not in df.columns:
        for alt in ("seq_pos", "ref_pos"):
            if alt in df.columns:
                df = df.rename(columns={alt: "pos"})
                break

    return df


def decode_base(val):
    """Convert dtw.base value (letter or integer) to a single uppercase nucleotide."""
    if pd.isna(val):
        return "N"
    if isinstance(val, (int, float)):
        return _INT_TO_BASE.get(int(val), "N")
    s = str(val).strip().upper()
    return s[0] if s and s[0] in "ACGTN" else "N"


def build_base_lookup(df):
    """
    Build a {(ref, pos_0based): nucleotide} dict from dtw.base column.
    Returns None if the column is absent or has no usable values.
    """
    if "dtw_base" not in df.columns:
        return None
    sub = df[["ref", "pos", "dtw_base"]].dropna(subset=["dtw_base"])
    if sub.empty:
        return None
    lookup = {}
    for _, row in sub.iterrows():
        key = (row["ref"], int(row["pos"]))
        if key not in lookup:
            lookup[key] = decode_base(row["dtw_base"])
    print(f"dtw.base lookup built: {len(lookup)} unique positions")
    return lookup


def main():
    args = parse_args()

    print("=== Per-base Signal Deviation ===")
    print(f"Input TSV: {args.tsv}")
    print(f"Genome:    {args.genome}")
    print(f"Output:    {args.output}")
    print()

    # na_values=["*"]: uncalled4 marks DTW failures with "*" (js4007 known issue).
    # on_bad_lines="warn": rarely, two TSV lines are concatenated without a newline
    #   producing a row with too many fields — skip rather than abort (js4007 known issue).
    # EmptyDataError: uncalled4 writes a 0-byte file when no reads pass strand filtering.
    try:
        df = pd.read_csv(args.tsv, sep="\t", na_values=["*", "NA", "nan"],
                         on_bad_lines="warn")
    except pd.errors.EmptyDataError:
        df = pd.DataFrame()
    if len(df) == 0:
        print("Warning: empty TSV input — writing empty output file")
        with gzip.open(args.output, "wt"):
            pass
        return

    df = normalize_columns(df)

    if "ref" not in df.columns or "pos" not in df.columns:
        sys.exit(
            f"Error: could not find ref/pos columns in TSV. "
            f"Columns present: {list(df.columns)}"
        )
    if "dtw_model_diff" not in df.columns:
        sys.exit(
            f"Error: 'dtw.model_diff' column not found in TSV. "
            f"Columns present: {list(df.columns)}"
        )

    print(f"Rows loaded: {len(df)}")
    print(f"Chromosomes: {sorted(df['ref'].unique())}")

    # Try to build nucleotide lookup from dtw.base (avoids pysam FASTA dependency).
    base_lookup = build_base_lookup(df)
    fasta = None
    if base_lookup is None:
        print("dtw.base column absent or empty — falling back to pysam FASTA lookup")
        try:
            import pysam
            fasta = pysam.FastaFile(args.genome)
        except ImportError:
            sys.exit(
                "Error: dtw.base not available and pysam not installed. "
                "Install with: conda install -c bioconda pysam"
            )

    # Group by (chromosome, position) — positions in TSV are 0-based per uncalled4 convention
    grouped = (
        df.groupby(["ref", "pos"])["dtw_model_diff"]
        .agg(
            mean_dev="mean",
            coverage="count",
            q25=lambda x: x.quantile(0.25),
            q75=lambda x: x.quantile(0.75),
            q025=lambda x: x.quantile(0.025),
            q975=lambda x: x.quantile(0.975),
        )
        .reset_index()
    )

    written = 0
    with gzip.open(args.output, "wt") as out:
        for _, row in grouped.sort_values(["ref", "pos"]).iterrows():
            chrom = row["ref"]
            pos_0 = int(row["pos"])   # 0-based from uncalled4
            cov = int(row["coverage"])
            mean_dev = float(row["mean_dev"])
            q25 = float(row["q25"])
            q75 = float(row["q75"])
            q025 = float(row["q025"])
            q975 = float(row["q975"])

            if cov < args.min_cov:
                continue

            if base_lookup is not None:
                nt = base_lookup.get((chrom, pos_0), "N")
            else:
                try:
                    nt = fasta.fetch(chrom, pos_0, pos_0 + 1).upper()
                except (ValueError, KeyError):
                    nt = "N"
                if not nt:
                    nt = "N"

            pos_1 = pos_0 + 1  # convert to 1-based for output
            out.write(
                f"{chrom}\t{pos_1}\t{nt}\t{cov}\t{mean_dev:.6f}"
                f"\t{q25:.6f}\t{q75:.6f}\t{q025:.6f}\t{q975:.6f}\n"
            )
            written += 1

    if fasta is not None:
        fasta.close()
    print(f"Done. Positions written: {written}")


if __name__ == "__main__":
    main()
