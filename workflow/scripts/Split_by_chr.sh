#!/usr/bin/env python3

"""Split by Chromosome: Split a tab-delimited file by chromosome.

Splits input by the first column (chromosome) into individual gzipped files.
Uses Python file I/O throughout — macOS BSD awk silently drops output files
during parallel Snakemake execution, even with direct file redirection.

Usage:
    Split_by_chr.sh -i <input_file> [-d <output_dir>]

Required arguments:
    -i    Input file (tab-delimited with chromosome in column 1)
          Supported formats: .txt, .txt.gz, .bed, .bed.gz, .bg, .bg.gz

Optional arguments:
    -d    Output directory (default: same directory as input)
    -h    Show this help message

Output:
    Creates <prefix>_<chr>.<ext>.gz for each chromosome.

Example:
    Split_by_chr.sh -i data/perbase_error/sample_for.txt.gz
    # Creates: data/perbase_error/sample_for_chr1.txt.gz, ...
"""

import argparse
import gzip
import glob
import os
import sys


def main():
    parser = argparse.ArgumentParser(description="Split a tab-delimited file by chromosome")
    parser.add_argument("-i", required=True, help="Input file")
    parser.add_argument("-d", default=None, help="Output directory")
    args = parser.parse_args()

    input_path = args.i

    if not os.path.isfile(input_path):
        print(f"Error: Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    # Determine output prefix and extension
    basename = os.path.basename(input_path)
    compressed = basename.endswith(".gz")
    if compressed:
        basename = basename[:-3]

    ext = None
    for candidate in ("txt", "bed", "bg"):
        if basename.endswith(f".{candidate}"):
            ext = candidate
            prefix = basename[: -(len(candidate) + 1)]
            break

    if ext is None:
        print("Error: Input file must have .txt, .bed, or .bg extension (optionally gzipped)", file=sys.stderr)
        sys.exit(1)

    if args.d is not None:
        outdir = args.d
        os.makedirs(outdir, exist_ok=True)
    else:
        outdir = os.path.dirname(input_path)

    out_prefix = os.path.join(outdir, f"{prefix}_")

    print("=== Split by Chromosome ===")
    print(f"Input: {input_path}")
    print(f"Output prefix: {out_prefix}")
    print(f"Output extension: {ext}")
    print()

    # Clean stale files from previous runs
    print("Cleaning output directory...")
    for stale in glob.glob(f"{out_prefix}*.{ext}") + glob.glob(f"{out_prefix}*.{ext}.gz"):
        os.remove(stale)

    # Split file by chromosome, writing gzipped output directly
    print("Splitting file by chromosome...")
    opener = gzip.open if compressed else open
    files = {}

    with opener(input_path, "rt") as fh:
        for line in fh:
            chr_name = line.split("\t", 1)[0]
            if chr_name not in files:
                out_path = f"{out_prefix}{chr_name}.{ext}.gz"
                files[chr_name] = gzip.open(out_path, "wt")
            files[chr_name].write(line)

    for f in files.values():
        f.flush()
        f.close()

    print(f"Split into {len(files)} chromosome files")

    # Verify every chromosome has its output file
    print("Verifying output files...")
    missing = []
    for chr_name in sorted(files):
        out_path = f"{out_prefix}{chr_name}.{ext}.gz"
        if not os.path.isfile(out_path):
            missing.append(out_path)
            print(f"  MISSING: {out_path}", file=sys.stderr)

    if missing:
        print(f"ERROR: {len(missing)} of {len(files)} chromosome files missing!", file=sys.stderr)
        sys.exit(1)

    print(f"All {len(files)} chromosome files verified.")
    print("Done.")


if __name__ == "__main__":
    main()
