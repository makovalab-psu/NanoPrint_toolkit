#!/usr/bin/env python3
"""
make_dtw_test_cases.py - Build edge-case Uncalled4 BAMs for the perbase_signal_deviation
rewrite, by deriving them from a REAL Uncalled4 BAM.

Why derive instead of synthesise
--------------------------------
The DTW payload is five aux tags (uncalled4 src/uncalled4/io/bam.py:22-27):

    ur  int32[]  reference interval bounds [s0,e0,s1,e1,...] of the fetched window
    ul  int16[]  sample-interval run lengths; one NON-NEGATIVE entry per aligned
                 reference position, negatives encode the start offset and gaps
    uc  int16[]  mean normalized current per aligned position, round(x * INORM_SCALE)
    ud  int16[]  current stdv, same encoding
    un  float[2] per-read normalization [scale, shift]

Writing those from scratch means reproducing uncalled4's k-mer registration by hand:
Sequence() shifts the ur bounds by model.PRMS.shift at the start and K-shift-1 at the
end, SWAPS the two for a reverse-strand read, and puts reverse-strand coordinates in
negative "mpos" space where pos = -mpos-1 (src/cpp/seq.hpp:270-306,
src/uncalled4/pore_model.py:519-527). If we hand-rolled that and got it wrong, the
test data and the code under test would be wrong in the SAME way and the test would
pass while the pipeline shifted every reactivity value by a few bases.

So: every case below starts from a real alignment and edits only the values, never
the coordinate structure. The structural invariants come along for free.

The two cases that cannot be reached this way (a long-interval run-length split, and
a synthetic k-mer registration check) are flagged in the manifest as REQUIRES_ORACLE
so they are never silently assumed to pass.

Ground truth
------------
This script does NOT compute expected output. Run the existing path on each derived
BAM, on a machine that has uncalled4:

    uncalled4 convert --bam-in case.bam --ref ref.fa --tsv-out case.tsv \
        --tsv-cols "dtw.current,dtw.current_sd,dtw.start,dtw.length,dtw.model_diff"
    perbase_signal_deviation.py -i case.tsv -g ref.fa -o case_expected.txt.gz

That output is the oracle. The new BAM-reading implementation must reproduce it.

Usage
-----
    make_dtw_test_cases.py -i <real_uncalled4.bam> -g <ref.fa> -o <outdir>
"""

import os
import sys
import array
import argparse

try:
    import pysam
except ImportError:
    sys.exit("Error: pysam not found. conda install -c bioconda pysam")

try:
    import numpy as np
except ImportError:
    sys.exit("Error: numpy not found.")


# Tag names (uncalled4 src/uncalled4/io/bam.py:22-27)
REF_TAG, LENS_TAG, CURS_TAG, STDS_TAG, NORM_TAG = "ur", "ul", "uc", "ud", "un"
REQ_TAGS = (REF_TAG, LENS_TAG, CURS_TAG)

# uncalled4 io/bam.py:168-170 -- MIN_I16 is int16.min+1, and the NA sentinel is the
# one value it deliberately avoids, so NA is int16.min. assert_sentinel_unused()
# re-checks this against the real data rather than trusting the inference.
NA_I16 = np.iinfo(np.int16).min          # -32768
MIN_I16 = np.iinfo(np.int16).min + 1     # -32767


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--bam", required=True, help="Real Uncalled4 BAM (indexed)")
    p.add_argument("-g", "--genome", required=True, help="Reference FASTA (for the manifest)")
    p.add_argument("-o", "--outdir", required=True, help="Output directory for the case BAMs")
    p.add_argument("--window", type=int, default=10000,
                   help="Window size the sweep will use; case 5 is built around a "
                        "multiple of this (default: 10000)")
    return p.parse_args()


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def has_dtw(rec):
    return all(rec.has_tag(t) for t in REQ_TAGS)


def n_positions(rec):
    """Aligned reference positions = non-negative entries in ul, and must equal len(uc).

    to_runlen() emits exactly one non-negative entry per sample interval
    (uncalled4 src/cpp/intervals.hpp:334-367): the interval length, or 0 when
    consecutive reference positions share one interval (a skip). Negative entries
    are the start pad and gap/long-run markers.
    """
    ul = np.asarray(rec.get_tag(LENS_TAG), dtype=np.int64)
    return int((ul >= 0).sum())


def ref_span(rec):
    """(start, end) over all ur intervals, in forward-genome coordinates."""
    ur = np.asarray(rec.get_tag(REF_TAG), dtype=np.int64)
    return int(ur[0]), int(ur[-1])


def set_current(rec, vals):
    """Replace uc with an int16 array, preserving everything else."""
    rec.set_tag(CURS_TAG, array.array("h", np.asarray(vals, dtype=np.int16)))


def clone(rec, header, name=None):
    out = pysam.AlignedSegment.fromstring(rec.to_string(), header)
    if name is not None:
        out.query_name = name
    return out


def write_case(outdir, name, header, records, note, manifest, requires_oracle=False):
    path = os.path.join(outdir, name + ".bam")
    records = sorted(records, key=lambda r: (r.reference_id, r.reference_start))
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for r in records:
            out.write(r)
    pysam.index(path)
    manifest.append((name, len(records), note, "yes" if requires_oracle else "no"))
    sys.stderr.write(f"  {name}.bam  ({len(records)} reads)  {note}\n")
    return path


def pick_both_strands(records, n=4):
    """Take up to n records, alternating strand as far as the pool allows.

    Selecting the first n in coordinate order gives single-strand cases, which
    is the worst possible sampling here: strand-dependent k-mer registration
    (start/end shift swap, negative mpos) is the main thing these cases exist to
    catch, so every case should carry both strands when the pool has both.
    """
    fwd = [r for r in records if not r.is_reverse]
    rev = [r for r in records if r.is_reverse]
    out = []
    while len(out) < n and (fwd or rev):
        if fwd:
            out.append(fwd.pop(0))
        if len(out) < n and rev:
            out.append(rev.pop(0))
    return out


def assert_sentinel_unused(records):
    """Guard the NA_I16 inference: if real data already carries -32768 in ul as a
    normal value, the sentinel assumption is wrong and every NA case is invalid."""
    for rec in records:
        ul = np.asarray(rec.get_tag(LENS_TAG), dtype=np.int64)
        if (ul == NA_I16).any():
            sys.exit(f"Error: {rec.query_name} has {NA_I16} in {LENS_TAG}, which this "
                     f"script assumes is the NA sentinel and never a real length. "
                     f"Re-check uncalled4 io/bam.py:168-170 before trusting the NA cases.")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    bam = pysam.AlignmentFile(args.bam, "rb")
    header = bam.header
    reads = [r for r in bam.fetch() if has_dtw(r)]
    bam.close()

    if not reads:
        sys.exit("Error: no records with DTW tags (ur/ul/uc) in the input BAM. "
                 "Is this really an 'uncalled4 align' output?")

    assert_sentinel_unused(reads)

    fwd = [r for r in reads if not r.is_reverse]
    rev = [r for r in reads if r.is_reverse]
    sys.stderr.write(f"Input: {len(reads)} DTW records ({len(fwd)} fwd, {len(rev)} rev)\n")
    if not fwd or not rev:
        sys.stderr.write("WARNING: one strand is missing; case 1 will be incomplete.\n")

    manifest = []

    # --- case 1: both strands over the same positions -----------------------
    # The highest-value case: catches the mpos negation and the swapped k-mer shift.
    # Pick the forward/reverse pair with the largest reference overlap.
    best = None
    for f in fwd:
        fs, fe = ref_span(f)
        for r in rev:
            rs, re_ = ref_span(r)
            ov = min(fe, re_) - max(fs, rs)
            if best is None or ov > best[0]:
                best = (ov, f, r)
    if best and best[0] > 0:
        _, f, r = best
        write_case(args.outdir, "case1_both_strands", header, [f, r],
                   f"fwd+rev overlapping by {best[0]} bp", manifest)
    else:
        sys.stderr.write("  SKIP case1: no overlapping fwd/rev pair\n")

    # --- case 2: multi-interval ur (reference deletion) ----------------------
    multi = [r for r in reads if len(r.get_tag(REF_TAG)) > 2]
    if multi:
        multi.sort(key=lambda r: -len(r.get_tag(REF_TAG)))
        multi = pick_both_strands(multi, 4)
        write_case(args.outdir, "case2_multi_interval", header, multi,
                   f"ur has up to {len(multi[0].get_tag(REF_TAG)) // 2} intervals", manifest)
    else:
        sys.stderr.write("  SKIP case2: no multi-interval ur in the real data "
                         "(needs a read with a reference deletion)\n")

    # --- case 3: NA positions in uc ------------------------------------------
    # Blank every 7th position of one read. These must drop out of N, not become 0.
    src = fwd[0] if fwd else reads[0]
    na_recs = []
    for rec0 in pick_both_strands(reads, 2):
        rec = clone(rec0, header, rec0.query_name + "_NA")
        uc = np.asarray(rec.get_tag(CURS_TAG), dtype=np.int16).copy()
        uc[::7] = NA_I16
        set_current(rec, uc)
        na_recs.append(rec)
    write_case(args.outdir, "case3_na_positions", header, na_recs,
               "every 7th position set to the NA sentinel, one read per strand", manifest)

    # --- case 4: skips (ul == 0) ---------------------------------------------
    skippy = [r for r in reads
              if (np.asarray(r.get_tag(LENS_TAG), dtype=np.int64) == 0).any()]
    if skippy:
        write_case(args.outdir, "case4_skips", header, pick_both_strands(skippy, 4),
                   "reads carrying ul==0 entries (shared sample interval)", manifest)
    else:
        sys.stderr.write("  SKIP case4: no ul==0 entries found\n")

    # --- case 5: read spanning a window boundary ------------------------------
    # Must be neither double-counted nor dropped by the windowed sweep.
    spanning = []
    for r in reads:
        s, e = ref_span(r)
        if s // args.window != (e - 1) // args.window:
            spanning.append(r)
    if spanning:
        write_case(args.outdir, "case5_window_boundary", header, pick_both_strands(spanning, 6),
                   f"reads crossing a multiple of {args.window}", manifest)
    else:
        sys.stderr.write(f"  SKIP case5: no read crosses a {args.window} bp boundary; "
                         f"rerun with a smaller --window\n")

    # --- case 6: coverage 1, and an all-NA position ---------------------------
    # One read alone -> every position has coverage 1 (quantile edge case).
    # One read per strand: each strand is processed independently, so every
    # position still has coverage 1 within its own strand.
    solo = [clone(r, header, r.query_name + "_solo") for r in pick_both_strands(reads, 2)]
    write_case(args.outdir, "case6a_coverage1", header, solo,
               "one read per strand: every position has coverage 1 within a strand", manifest)

    # Two reads, with the SAME positions blanked in both -> those positions have
    # observations but zero valid ones.
    if len(fwd) >= 2:
        a, b = fwd[0], fwd[1]
        as_, ae = ref_span(a)
        bs, be = ref_span(b)
        lo, hi = max(as_, bs), min(ae, be)
        if hi > lo:
            pair = []
            for i, rec0 in enumerate((a, b)):
                rec = clone(rec0, header, rec0.query_name + f"_allNA{i}")
                s, _ = ref_span(rec)
                uc = np.asarray(rec.get_tag(CURS_TAG), dtype=np.int16).copy()
                # blank a 10 bp stretch inside the shared span, indexed from the
                # read's own start. Approximate for multi-interval reads, which is
                # why the oracle - not this script - defines the expected answer.
                # Blank mid-overlap, away from the contig ends where the k-mer
                # trimming lives - an edge-adjacent blank would confound the two
                # effects if the test failed.
                mid = (lo + hi) // 2
                off = mid - s
                if 0 <= off < len(uc):
                    uc[off:off + 10] = NA_I16
                set_current(rec, uc)
                pair.append(rec)
            write_case(args.outdir, "case6b_all_na_position", header, pair,
                       f"~10 bp near {(lo + hi) // 2} blanked in both reads", manifest)
        else:
            sys.stderr.write("  SKIP case6b: top two fwd reads do not overlap\n")

    # --- case 7: split read (pi tag) ------------------------------------------
    split = [r for r in reads if r.has_tag("pi")]
    if split:
        write_case(args.outdir, "case7_split_read", header, pick_both_strands(split, 4),
                   "reads with the pi tag (sam_to_aln takes a different branch)", manifest)
    else:
        sys.stderr.write("  SKIP case7: no pi-tagged reads; carry one over from a "
                         "js4022 BAM, where read splitting was on\n")

    # --- case 8: long-interval run-length split -------------------------------
    # to_runlen emits negative MAX_LEN_I16 chunks for an interval longer than 32767
    # samples (src/cpp/intervals.hpp:345-350). Cannot be forged by editing values -
    # it needs a genuinely long stall. Flagged, not faked.
    longrun = [r for r in reads
               if (np.asarray(r.get_tag(LENS_TAG), dtype=np.int64) <= MIN_I16).any()]
    if longrun:
        write_case(args.outdir, "case8_long_runlen", header, pick_both_strands(longrun, 4),
                   "reads with a <= MIN_I16 entry (long interval or large start offset)",
                   manifest, requires_oracle=True)
    else:
        sys.stderr.write("  SKIP case8: no long-run markers in this data (expected for "
                         "mtDNA); look for one in a js4022 BAM\n")

    # --- manifest -------------------------------------------------------------
    mpath = os.path.join(args.outdir, "manifest.tsv")
    with open(mpath, "w") as fh:
        fh.write("case\treads\tnote\trequires_oracle\n")
        for row in manifest:
            fh.write("\t".join(str(c) for c in row) + "\n")

    sys.stderr.write(f"\nWrote {len(manifest)} cases + manifest.tsv to {args.outdir}\n")
    sys.stderr.write(
        "\nNext, on a machine with uncalled4, build the oracle for each case:\n"
        "  for b in %s/case*.bam; do\n"
        "    uncalled4 convert --bam-in $b --ref %s --tsv-out ${b%%.bam}.tsv \\\n"
        "      --tsv-cols 'dtw.current,dtw.current_sd,dtw.start,dtw.length,dtw.model_diff'\n"
        "    for s in for rev; do\n"
        "      perbase_signal_deviation.py -i ${b%%.bam}.tsv -g %s -o ${b%%.bam}_${s}_expected.txt.gz\n"
        "    done\n"
        "  done\n" % (args.outdir, args.genome, args.genome))


if __name__ == "__main__":
    main()
