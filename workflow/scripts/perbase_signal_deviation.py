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
  - dtw_model_diff = observed current - model current, in NORMALIZED units, not pA.
    uncalled4 computes `np.array(self.current) - self.seq.current` (tracks.py:213), i.e.
    observed MINUS predicted, so positive = observed current HIGHER than the pore model.
    Both terms live in the normalized space the signal was scaled into before DTW, so the
    difference is unitless; multiply by the model's pa_stdv (~23.1 pA for R10.4.1) for pA.

MEMORY: the TSV holds one row per aligned base per read, so at WGS depth it is far
larger than RAM — js4022 is ~9.5e9 rows (~570 GB) for a single sample strand at 2000x
over a 4.2 Mb genome. Reading it whole, as this script used to with pandas.read_csv,
is impossible there. It now works in two passes over bounded memory:

  Pass 1  Stream the TSV line by line and bin every observation into a fixed-width
          genome window (default 10 kb) under a temp directory. Rows are buffered in
          compact array.array buffers and flushed in batches, so the resident set is
          the buffer, not the file.
  Pass 2  Load one window at a time and compute the per-position statistics for it
          with numpy, writing results in coordinate order.

Peak memory is therefore set by the busiest single window — coverage x window size x
12 bytes, about 120 MB for a 10 kb window at 1000x per strand — plus the pass 1
buffer. Lower --window for deeper data or a tighter memory limit; it changes nothing
but the peak, since windows are cut on genome coordinates and every observation at a
position lands in the same window.

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
    5. Mean signal deviation (mean dtw.model_diff, normalized units)
    6. Q25 — 0.25 quantile of dtw.model_diff (lower 50% CI bound, normalized)
    7. Q75 — 0.75 quantile of dtw.model_diff (upper 50% CI bound, normalized)
    8. Q025 — 0.025 quantile of dtw.model_diff (lower 95% CI bound, normalized)
    9. Q975 — 0.975 quantile of dtw.model_diff (upper 95% CI bound, normalized)
   10. Mean squared deviation — mean(dtw.model_diff^2) = sum(dtw.model_diff^2) / N (normalized^2)

Quantiles use linear interpolation between order statistics, the default of both
numpy.percentile and pandas.Series.quantile, so output matches the previous
pandas implementation.

Usage:
    perbase_signal_deviation.py -i <dtw.tsv> -g <genome.fa> -o <output.txt.gz>
      [--window 10000] [--buffer-rows 20000000] [--tmp-dir DIR] [-c MIN_COV]
"""

import os
import sys
import gzip
import shutil
import argparse
import tempfile
from array import array

try:
    import numpy as np
except ImportError:
    sys.exit("Error: numpy not found. Install with: conda install numpy")


# Integer encoding used by uncalled4's "binarized" dtw.base column (if applicable)
_INT_TO_BASE = {0: "A", 1: "C", 2: "G", 3: "T"}

# Values uncalled4 writes for a position where DTW failed, plus the usual null spellings.
_NULL_VALUES = {"*", "", "NA", "nan", "NaN", "None"}

# Record layout of the per-window temp files: position, value, and (only when the TSV
# carries a usable dtw.base column) the encoded nucleotide.
_POS_T = "i"      # int32  — genome coordinate, 0-based as uncalled4 writes it
_VAL_T = "d"      # float64 — keeps sums bit-identical to the old pandas version
_BASE_T = "B"     # uint8  — ord() of the nucleotide letter


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
    p.add_argument("-w", "--window", type=int, default=10000,
                   help="Genome window size in nt for binning (default: 10000). "
                        "Sets peak memory: one window's observations are held at once.")
    p.add_argument("--buffer-rows", type=int, default=20_000_000,
                   help="Rows buffered in memory before flushing to window files "
                        "(default: 20000000, about 240 MB)")
    p.add_argument("--tmp-dir", default=None,
                   help="Directory for the per-window temp files (default: alongside "
                        "the output). Needs room for roughly the size of the input TSV.")
    return p.parse_args()


def normalize_name(name):
    """Uncalled4 column name -> canonical name used here."""
    name = name.strip().replace(".", "_")
    if name in ("seq_name", "chr"):
        return "ref"
    if name in ("seq_pos", "ref_pos"):
        return "pos"
    return name


def decode_base(token):
    """Convert a dtw.base field (letter or integer) to a single uppercase nucleotide."""
    token = token.strip().upper()
    if not token or token in _NULL_VALUES:
        return "N"
    if token[0] in "ACGTN":
        return token[0]
    try:
        return _INT_TO_BASE.get(int(float(token)), "N")
    except ValueError:
        return "N"


class WindowWriter:
    """Bins observations into per-window temp files, buffering in compact arrays.

    One file per (chromosome, window). Buffers are array.array, so a buffered row
    costs 12 bytes rather than the ~60 a Python tuple of int and float would.
    """

    def __init__(self, tmp_dir, window, max_rows, with_bases):
        self.tmp_dir = tmp_dir
        self.window = window
        self.max_rows = max_rows
        self.with_bases = with_bases
        self.buffers = {}          # (ref, win) -> [pos array, val array, base array|None]
        self.buffered = 0
        self.paths = {}            # (ref, win) -> temp file path
        self.ref_names = {}        # ref -> filename-safe token
        self.flushes = 0

    def _path(self, key):
        path = self.paths.get(key)
        if path is None:
            ref, win = key
            token = self.ref_names.get(ref)
            if token is None:
                token = "ref%d" % len(self.ref_names)
                self.ref_names[ref] = token
            path = os.path.join(self.tmp_dir, "%s_%09d.bin" % (token, win))
            self.paths[key] = path
        return path

    def add(self, ref, pos, value, base_code):
        key = (ref, pos // self.window)
        buf = self.buffers.get(key)
        if buf is None:
            buf = [array(_POS_T), array(_VAL_T),
                   array(_BASE_T) if self.with_bases else None]
            self.buffers[key] = buf
        buf[0].append(pos)
        buf[1].append(value)
        if self.with_bases:
            buf[2].append(base_code)
        self.buffered += 1
        if self.buffered >= self.max_rows:
            self.flush()

    def flush(self):
        if not self.buffered:
            return
        for key, buf in self.buffers.items():
            if not buf[0]:
                continue
            with open(self._path(key), "ab") as fh:
                buf[0].tofile(fh)
                buf[1].tofile(fh)
                if self.with_bases:
                    buf[2].tofile(fh)
            # A window file is a sequence of independently written chunks; pass 2
            # reads each chunk back by its own row count, so chunks never mix.
            with open(self._path(key) + ".idx", "a") as idx:
                idx.write("%d\n" % len(buf[0]))
            del buf[0][:]
            del buf[1][:]
            if self.with_bases:
                del buf[2][:]
        self.buffers.clear()
        self.buffered = 0
        self.flushes += 1

    def windows(self):
        """(ref, win, path) for every window written, in coordinate order."""
        return sorted(self.paths.items(), key=lambda kv: (kv[0][0], kv[0][1]))


def read_window(path, with_bases):
    """Read one window file back into numpy arrays."""
    counts = [int(line) for line in open(path + ".idx")]
    pos_parts, val_parts, base_parts = [], [], []
    with open(path, "rb") as fh:
        for n in counts:
            pos_parts.append(np.fromfile(fh, dtype="<i4", count=n))
            val_parts.append(np.fromfile(fh, dtype="<f8", count=n))
            if with_bases:
                base_parts.append(np.fromfile(fh, dtype="u1", count=n))
    pos = np.concatenate(pos_parts) if len(pos_parts) > 1 else pos_parts[0]
    val = np.concatenate(val_parts) if len(val_parts) > 1 else val_parts[0]
    base = None
    if with_bases:
        base = np.concatenate(base_parts) if len(base_parts) > 1 else base_parts[0]
    return pos, val, base


def group_quantile(values_sorted, starts, counts, q):
    """Linear-interpolation quantile for each group of a value-sorted array.

    Groups are contiguous runs given by starts/counts, each sorted ascending. This is
    the same estimator as numpy.percentile and pandas.Series.quantile with their
    default method: h = (n - 1) * q, interpolating between order statistics h and h+1.

    The interpolation mirrors numpy's own `_lerp`, which switches to `b - (b - a) *
    (1 - t)` once t >= 0.5 for numerical stability. Using the plain `a + t * (b - a)`
    everywhere agrees to about 1e-6 — enough to move the last printed digit of a
    quantile column and make outputs from the two implementations diff.
    """
    h = (counts - 1) * q
    lo = np.floor(h).astype(np.int64)
    frac = h - lo
    lo_idx = starts + lo
    hi_idx = np.minimum(lo_idx + 1, starts + counts - 1)
    low = values_sorted[lo_idx]
    high = values_sorted[hi_idx]
    diff = high - low
    result = low + diff * frac
    upper = frac >= 0.5
    result[upper] = high[upper] - diff[upper] * (1 - frac[upper])
    return result


def main():
    args = parse_args()

    print("=== Per-base Signal Deviation ===")
    print(f"Input TSV: {args.tsv}")
    print(f"Genome:    {args.genome}")
    print(f"Output:    {args.output}")
    print(f"Window:    {args.window} nt")
    print()

    out_dir = os.path.dirname(os.path.abspath(args.output))
    os.makedirs(out_dir, exist_ok=True)
    tmp_parent = args.tmp_dir or out_dir
    os.makedirs(tmp_parent, exist_ok=True)
    tmp_dir = tempfile.mkdtemp(prefix=".perbase_signal_", dir=tmp_parent)

    try:
        written = run(args, tmp_dir)
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    print(f"Done. Positions written: {written}")


def run(args, tmp_dir):
    # ------------------------------------------------------------------
    # Pass 1: stream the TSV, binning observations into window files
    # ------------------------------------------------------------------
    opener = gzip.open if args.tsv.endswith(".gz") else open
    with opener(args.tsv, "rt") as fh:
        header = fh.readline()
        if not header.strip():
            print("Warning: empty TSV input — writing empty output file")
            with gzip.open(args.output, "wt"):
                pass
            return 0

        columns = [normalize_name(c) for c in header.rstrip("\n").split("\t")]
        index = {name: i for i, name in enumerate(columns)}

        if "ref" not in index or "pos" not in index:
            sys.exit(
                f"Error: could not find ref/pos columns in TSV. "
                f"Columns present: {columns}"
            )
        if "dtw_model_diff" not in index:
            sys.exit(
                f"Error: 'dtw.model_diff' column not found in TSV. "
                f"Columns present: {columns}"
            )

        i_ref = index["ref"]
        i_pos = index["pos"]
        i_val = index["dtw_model_diff"]
        i_base = index.get("dtw_base")
        n_cols = len(columns)
        with_bases = i_base is not None

        writer = WindowWriter(tmp_dir, args.window, args.buffer_rows, with_bases)
        rows = skipped_null = skipped_bad = 0
        base_seen = False

        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != n_cols:
                # Rarely uncalled4 writes two rows without the newline between them
                # (js4007 known issue). Skip and report, as the old on_bad_lines="warn".
                skipped_bad += 1
                continue
            raw = fields[i_val]
            if raw in _NULL_VALUES:
                skipped_null += 1
                continue
            try:
                value = float(raw)
                pos = int(fields[i_pos])
            except ValueError:
                skipped_bad += 1
                continue
            if value != value:      # NaN
                skipped_null += 1
                continue
            code = 0
            if with_bases:
                nt = decode_base(fields[i_base])
                if nt != "N":
                    base_seen = True
                code = ord(nt)
            writer.add(fields[i_ref], pos, value, code)
            rows += 1

        writer.flush()

    print(f"Rows binned: {rows}")
    if skipped_null:
        print(f"Rows skipped (no DTW value): {skipped_null}")
    if skipped_bad:
        print(f"Warning: rows skipped (malformed): {skipped_bad}")
    windows = writer.windows()
    print(f"Windows on disk: {len(windows)} (flushes: {writer.flushes})")

    if rows == 0:
        print("Warning: no usable rows — writing empty output file")
        with gzip.open(args.output, "wt"):
            pass
        return 0

    # dtw.base is only usable if the column existed AND held real nucleotides.
    with_bases = with_bases and base_seen
    fasta = None
    if not with_bases:
        print("dtw.base column absent or empty — falling back to pysam FASTA lookup")
        try:
            import pysam
            fasta = pysam.FastaFile(args.genome)
        except ImportError as exc:
            # Report the real exception. "import pysam" raises ImportError both when
            # pysam is absent and when it is present but its C extension cannot load
            # a shared library — the second case looks identical to the first, and
            # collapsing them into "not installed" sends you chasing the wrong bug
            # (notably when a submit node imports it fine and a compute node does not).
            sys.exit(
                f"Error: could not import pysam or open the reference.\n"
                f"  {type(exc).__name__}: {exc}\n"
                f"  python:  {sys.executable}\n"
                f"If the message names a missing .so, pysam is installed but its\n"
                f"C extension cannot load here — check the environment on the node\n"
                f"actually running the job, not the submit host.\n"
                f"If pysam is genuinely absent: pip install pysam"
            )

    # ------------------------------------------------------------------
    # Pass 2: one window at a time, in coordinate order
    # ------------------------------------------------------------------
    written = 0
    contig_cache = {}               # ref -> sequence of the current contig only

    with gzip.open(args.output, "wt") as out:
        for (ref, _win), path in windows:
            pos, val, base = read_window(path, writer.with_bases)

            # Sort by position, and by value within each position, so the per-group
            # statistics below are all index arithmetic on contiguous runs.
            order = np.lexsort((val, pos))
            pos = pos[order]
            val = val[order]
            if base is not None:
                base = base[order]

            boundaries = np.flatnonzero(np.diff(pos)) + 1
            starts = np.concatenate(([0], boundaries))
            counts = np.diff(np.concatenate((starts, [pos.size])))
            positions = pos[starts]

            sums = np.add.reduceat(val, starts)
            sums_sq = np.add.reduceat(val * val, starts)
            means = sums / counts
            mean_sqs = sums_sq / counts
            q25 = group_quantile(val, starts, counts, 0.25)
            q75 = group_quantile(val, starts, counts, 0.75)
            q025 = group_quantile(val, starts, counts, 0.025)
            q975 = group_quantile(val, starts, counts, 0.975)

            if base is not None:
                nts = base[starts]

            if ref not in contig_cache:
                contig_cache.clear()    # one contig resident at a time
                if fasta is not None:
                    try:
                        contig_cache[ref] = fasta.fetch(ref).upper()
                    except (ValueError, KeyError):
                        contig_cache[ref] = ""
                else:
                    contig_cache[ref] = ""
            sequence = contig_cache[ref]

            for k in range(positions.size):
                cov = int(counts[k])
                if cov < args.min_cov:
                    continue
                pos_0 = int(positions[k])   # 0-based from uncalled4
                if base is not None:
                    nt = chr(int(nts[k]))
                else:
                    nt = sequence[pos_0] if 0 <= pos_0 < len(sequence) else "N"
                    if nt not in "ACGTN":
                        nt = "N"
                out.write(
                    f"{ref}\t{pos_0 + 1}\t{nt}\t{cov}\t{means[k]:.6f}"
                    f"\t{q25[k]:.6f}\t{q75[k]:.6f}\t{q025[k]:.6f}\t{q975[k]:.6f}"
                    f"\t{mean_sqs[k]:.6f}\n"
                )
                written += 1

    if fasta is not None:
        fasta.close()
    return written


if __name__ == "__main__":
    main()
