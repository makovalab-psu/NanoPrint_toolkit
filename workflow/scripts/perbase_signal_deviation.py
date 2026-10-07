#!/usr/bin/env python3
"""
perbase_signal_deviation.py - Per-base pore model signal deviation, read straight
from an Uncalled4 signal-alignment BAM.

Reads the DTW tags out of the BAM with uncalled4's own decoder. The intermediate
`uncalled4 convert` TSV is no longer involved: it held one row per aligned base per
read (~570 GB for a single js4022 sample strand), and building it cost more than the
statistics it fed. `uncalled4_convert_tsv` is still a Snakemake rule for small
datasets and for producing the oracle these outputs are tested against - it is just
no longer on the path to the signal tracks.

Why uncalled4's decoder and not our own tag parsing
---------------------------------------------------
The five DTW tags (ur/ul/uc/ud/un) are cheap to read, but turning them into genome
coordinates is not: Sequence() shifts the ur bounds by model.PRMS.shift at the start
and k-shift-1 at the end, SWAPS those two for a reverse-strand read, and places
reverse-strand coordinates in negative "mpos" space where pos = -mpos-1
(uncalled4 src/cpp/seq.hpp:270-306, src/uncalled4/pore_model.py:519-527). Getting
that wrong shifts every reactivity value by a few bases while the track still looks
entirely plausible. So the coordinates and the model_diff values come from
uncalled4; this script only accumulates them.

Streaming, in one pass
----------------------
The input BAM is coordinate-sorted, so observations arrive in genome order. Each
read's (position, deviation) pairs are binned into fixed-width genome windows, and a
window is finalised as soon as no later read can reach it - reads arrive sorted by
start, so window w is done once a read starts at or past its end. Nothing is written
to disk, and peak memory is the windows currently open: window size x per-strand
coverage x 12 bytes, about 120 MB for a 10 kb window at 1000x. Lower --window for
deeper data; it changes the peak, never the result.

For comparison, the path this replaces wrote a ~570 GB TSV, a strand-filtered copy of
the BAM, and its own binary window files.

Output format (tab-delimited, gzipped) - unchanged, so old and new outputs diff:
    1. Chromosome
    2. Position (1-based)
    3. Nucleotide
    4. Coverage (reads contributing to this position)
    5. Mean signal deviation (mean dtw.model_diff, normalized units)
    6. Q25   - 0.25 quantile   (lower 50% CI bound)
    7. Q75   - 0.75 quantile   (upper 50% CI bound)
    8. Q025  - 0.025 quantile  (lower 95% CI bound)
    9. Q975  - 0.975 quantile  (upper 95% CI bound)
   10. Mean squared deviation - mean(dtw.model_diff^2), normalized^2

dtw.model_diff is observed minus pore-model current in NORMALIZED units, not pA
(uncalled4 tracks.py:212-213). Multiply by the model's pa_stdv - recorded in the
BAM's own "@CO UNC:" header, 23.11 for R10.4.1 - to get pA.

Usage:
    perbase_signal_deviation.py -i <uncalled4.bam> -g <genome.fa> -s <for|rev> \
        -o <output.txt.gz> [--window 10000] [-c MIN_COV]
"""

import os
import sys
import gzip
import argparse

try:
    import numpy as np
except ImportError:
    sys.exit("Error: numpy not found. Install with: conda install numpy")


def parse_args():
    p = argparse.ArgumentParser(
        description="Per-base signal deviation from an Uncalled4 BAM",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    p.add_argument("-i", "--bam", required=True,
                   help="Input Uncalled4 BAM (coordinate-sorted and indexed)")
    p.add_argument("-g", "--genome", required=True,
                   help="Reference genome FASTA. Passed to uncalled4 as --ref (it needs "
                        "the reference to look up the pore model k-mer behind "
                        "dtw.model_diff) and used for nucleotide identity.")
    p.add_argument("-s", "--strand", required=True, choices=["for", "rev"],
                   help="Strand to process. Replaces the samtools pre-filter the old "
                        "convert-based rule needed: reads are filtered in memory.")
    p.add_argument("-o", "--output", required=True, help="Output file (.txt.gz)")
    p.add_argument("-c", "--min-cov", type=int, default=1,
                   help="Minimum coverage to emit a position (default: 1)")
    p.add_argument("-w", "--window", type=int, default=10000,
                   help="Genome window size in nt (default: 10000). Sets peak memory: "
                        "one window's observations are held at once.")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Statistics - unchanged from the TSV implementation, so outputs stay comparable
# ---------------------------------------------------------------------------

def group_quantile(values_sorted, starts, counts, q):
    """Linear-interpolation quantile for each group of a value-sorted array.

    Groups are contiguous runs given by starts/counts, each sorted ascending. Same
    estimator as numpy.percentile and pandas.Series.quantile with their default
    method: h = (n - 1) * q, interpolating between order statistics h and h+1.

    The interpolation mirrors numpy's own `_lerp`, which switches to
    `b - (b - a) * (1 - t)` once t >= 0.5 for numerical stability. Using the plain
    `a + t * (b - a)` everywhere agrees to about 1e-6 - enough to move the last
    printed digit of a quantile column and make the two implementations diff.
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


def format_window(ref, pos, val, sequence, min_cov):
    """Per-position statistics for one window's observations -> output lines.

    pos/val are all observations that fell in this window, in any order.
    """
    if pos.size == 0:
        return []

    # Sort by position, and by value within each position, so every statistic below
    # is index arithmetic on contiguous runs.
    order = np.lexsort((val, pos))
    pos = pos[order]
    val = val[order]

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

    lines = []
    for k in range(positions.size):
        cov = int(counts[k])
        if cov < min_cov:
            continue
        pos_0 = int(positions[k])       # 0-based, as uncalled4 reports it
        nt = sequence[pos_0] if 0 <= pos_0 < len(sequence) else "N"
        if nt not in "ACGTN":
            nt = "N"
        lines.append(
            f"{ref}\t{pos_0 + 1}\t{nt}\t{cov}\t{means[k]:.6f}"
            f"\t{q25[k]:.6f}\t{q75[k]:.6f}\t{q025[k]:.6f}\t{q975[k]:.6f}"
            f"\t{mean_sqs[k]:.6f}\n"
        )
    return lines


# ---------------------------------------------------------------------------
# Streaming window accumulator
# ---------------------------------------------------------------------------

class WindowAccumulator:
    """Bins observations into genome windows and finalises them in coordinate order.

    Relies on the input being coordinate-sorted: once a read starts at or beyond the
    end of window w, nothing can be added to w again. Windows below the watermark are
    handed to `emit` and freed.

    Arriving below the watermark would mean silently dropping observations, so it is
    a hard error rather than a warning - a track that is quietly missing a slice of
    its data looks completely normal downstream.
    """

    def __init__(self, window, emit):
        self.window = window
        self.emit = emit                # emit(ref, win_idx, pos_array, val_array)
        self.ref = None
        self.bufs = {}                  # win_idx -> ([pos arrays], [val arrays])
        self.watermark = 0              # lowest window index still open
        self.max_open = 0

    def _flush(self, below=None):
        keys = sorted(self.bufs) if below is None else \
               sorted(k for k in self.bufs if k < below)
        for k in keys:
            pos_parts, val_parts = self.bufs.pop(k)
            pos = np.concatenate(pos_parts) if len(pos_parts) > 1 else pos_parts[0]
            val = np.concatenate(val_parts) if len(val_parts) > 1 else val_parts[0]
            self.emit(self.ref, k, pos, val)
        if below is not None:
            self.watermark = max(self.watermark, below)

    def set_ref(self, ref):
        """Contig change: everything open belongs to the old contig."""
        if ref != self.ref:
            self._flush()
            self.ref = ref
            self.watermark = 0

    def add(self, pos, val, read_start):
        """Add one read's observations, then close windows it has moved past."""
        if pos.size:
            wi = pos // self.window
            lowest = int(wi.min())
            if lowest < self.watermark:
                raise RuntimeError(
                    f"Observation at {self.ref}:{int(pos.min())} falls in window "
                    f"{lowest}, already finalised (watermark {self.watermark}). The "
                    f"input BAM is not coordinate-sorted, or a read reports positions "
                    f"before its own alignment start. Sort the BAM and rerun."
                )
            # np.unique over a read's windows: a read touches (span / window) + 1 of
            # them, so this loop is one or two iterations in practice.
            for w in np.unique(wi):
                m = wi == w
                buf = self.bufs.get(int(w))
                if buf is None:
                    buf = ([], [])
                    self.bufs[int(w)] = buf
                buf[0].append(pos[m])
                buf[1].append(val[m])
            self.max_open = max(self.max_open, len(self.bufs))

        self._flush(below=read_start // self.window)

    def close(self):
        self._flush()


# ---------------------------------------------------------------------------
# The only uncalled4-dependent part
# ---------------------------------------------------------------------------

def iter_read_observations(bam_path, genome, want_reverse):
    """Yield (ref_name, read_start, positions, deviations) per alignment.

    positions are 0-based genome coordinates, deviations are dtw.model_diff with the
    NaN entries (positions where DTW failed) already removed. Both come from
    uncalled4's own decoder, so the k-mer registration and the reverse-strand
    coordinate flip are its code, not ours.
    """
    try:
        from uncalled4 import Config, Tracks
    except ImportError as exc:
        sys.exit(
            f"Error: could not import uncalled4.\n"
            f"  {type(exc).__name__}: {exc}\n"
            f"  python:  {sys.executable}\n"
            f"Install with: pip install setuptools==69.5.1 && pip install uncalled4\n"
            f"If the message names a missing .so, uncalled4 is installed but its C\n"
            f"extension cannot load here - check the environment on the node actually\n"
            f"running the job, not the submit host."
        )

    conf = Config()
    conf.tracks.io.bam_in = [bam_path]
    conf.tracks.ref = genome
    # "dtw" only. Adding "moves" makes sam_to_aln want the basecaller move table and,
    # through it, the read index - this script never opens a pod5.
    conf.tracks.layers = ["dtw"]
    conf.read_index.load_signal = False

    tracks = Tracks(conf=conf)
    bam_in = tracks.bam_in
    if bam_in is None:
        sys.exit(f"Error: uncalled4 did not open {bam_path} as a BAM input.")

    n_reads = n_skipped = n_nodtw = 0
    try:
        for sam in bam_in.iter_sam():
            # Filter before decoding: the other strand costs nothing this way.
            if bool(sam.is_reverse) != want_reverse:
                continue
            aln = bam_in.sam_to_aln(sam, load_moves=False)
            if aln is None:
                n_nodtw += 1
                continue

            pos = np.asarray(aln.seq.pos, dtype=np.int64)
            val = np.asarray(aln.dtw.model_diff, dtype=np.float64)

            if pos.size != val.size:
                # Never observed, but a silent mismatch here would pair deviations
                # with the wrong coordinates, which is the exact failure this design
                # exists to avoid. Refuse rather than guess at the alignment.
                sys.exit(
                    f"Error: read {sam.query_name} has {pos.size} reference positions "
                    f"but {val.size} DTW values. Refusing to guess how they line up."
                )

            keep = ~np.isnan(val)
            n_skipped += int((~keep).sum())
            n_reads += 1
            yield sam.reference_name, sam.reference_start, pos[keep], val[keep]
    finally:
        tracks.close()

    sys.stderr.write(f"Reads decoded: {n_reads}\n")
    if n_skipped:
        sys.stderr.write(f"Positions skipped (no DTW value): {n_skipped}\n")
    if n_nodtw:
        sys.stderr.write(f"Records skipped (no DTW tags): {n_nodtw}\n")


# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    print("=== Per-base Signal Deviation (direct from Uncalled4 BAM) ===")
    print(f"Input BAM: {args.bam}")
    print(f"Genome:    {args.genome}")
    print(f"Strand:    {args.strand}")
    print(f"Output:    {args.output}")
    print(f"Window:    {args.window} nt")
    print()

    if not os.path.exists(args.bam):
        sys.exit(f"Error: input BAM not found: {args.bam}")

    out_dir = os.path.dirname(os.path.abspath(args.output))
    os.makedirs(out_dir, exist_ok=True)

    try:
        import pysam
        fasta = pysam.FastaFile(args.genome)
    except ImportError as exc:
        # Report the real exception. "import pysam" raises ImportError both when pysam
        # is absent and when its C extension cannot load a shared library; collapsing
        # the two into "not installed" sends you chasing the wrong bug, notably when a
        # submit node imports it fine and a compute node does not.
        sys.exit(
            f"Error: could not import pysam.\n"
            f"  {type(exc).__name__}: {exc}\n"
            f"  python:  {sys.executable}\n"
            f"If the message names a missing .so, pysam is installed but its C\n"
            f"extension cannot load here. If it is genuinely absent: pip install pysam"
        )

    want_reverse = (args.strand == "rev")
    written = 0
    contig_cache = {}               # one contig resident at a time

    with gzip.open(args.output, "wt") as out:

        def emit(ref, _win_idx, pos, val):
            nonlocal written
            if ref not in contig_cache:
                contig_cache.clear()
                try:
                    contig_cache[ref] = fasta.fetch(ref).upper()
                except (ValueError, KeyError):
                    sys.stderr.write(
                        f"Warning: {ref} not found in {args.genome}; "
                        f"nucleotides for it will be N\n")
                    contig_cache[ref] = ""
            for line in format_window(ref, pos, val, contig_cache[ref], args.min_cov):
                out.write(line)
                written += 1

        acc = WindowAccumulator(args.window, emit)
        for ref, read_start, pos, val in iter_read_observations(
                args.bam, args.genome, want_reverse):
            acc.set_ref(ref)
            acc.add(pos, val, read_start)
        acc.close()

    fasta.close()

    if written == 0:
        print("Warning: no positions written. This is expected only when the sample "
              "has no reads on this strand.")
    print(f"Done. Positions written: {written}")
    print(f"Peak windows held open: {acc.max_open}")


if __name__ == "__main__":
    main()
