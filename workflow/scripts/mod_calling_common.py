#!/usr/bin/env python3
"""
mod_calling_common.py - shared core for direct modification calling by sample
comparison (phases 2c/3c). Imported by kmer_signal_model.py and perbase_mod_calls.py.

The method is Rembo's (P2-seq preprint, Methods "Rembo pipeline"), with the test
itself taken from Tombo's model_sample_compare
(tombo_stats.py compute_sample_compare_read_stats):

  1. Canonical model, from the UNTREATED sample: pool the per-read signal level of
     every k-mer observation whose k-mer aligned with no mismatch or indel, drop
     observations more than --mad median absolute deviations from the k-mer median,
     keep n / mean / SD per k-mer.
  2. Test, on any sample: z = (level - mean_kmer) / SD_kmer, two-sided Gaussian
     p-value, one per read per reference position. Optionally combined along the
     read by Fisher's method over +/- lag positions (Tombo's --fishers-method-context).
  3. A base is "modified" at p < threshold; modified fraction = modified reads /
     scored reads at that position.

Two backends supply step 1 and 2 with the same thing - (position, level, k-mer) per
read - from different signal-to-reference mappings:

  uncalled4   dtw.current and seq.kmer decoded from the Uncalled4 BAM, exactly as
              perbase_signal_deviation.py decodes dtw.model_diff. No pod5.
  remora      Remora API on pod5 + BAM: move table (optionally refined against a
              k-mer level table), per-base trimmed mean. This is the route Rembo
              itself takes.

Levels from the two are in different units and registered to different bases of
the k-mer, so a model is only ever applied to data from the backend that built it;
load_model() callers check this.

K-mers are in READ orientation (5'->3' as the strand went through the pore) and
encoded base-4, A=0 C=1 G=2 T=3, 5' base most significant - the same encoding
uncalled4 uses (pore_model.hpp str_to_kmer), verified at start-up.
"""

import sys
import zlib
import gzip
import math
from statistics import NormalDist

import numpy as np

BACKENDS = ("uncalled4", "remora")

_BASE_LUT = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate("ACGT"):
    _BASE_LUT[ord(_b)] = _i
    _BASE_LUT[ord(_b.lower())] = _i

SMALLEST_P = 1e-300


# ---------------------------------------------------------------------------
# Sequence helpers
# ---------------------------------------------------------------------------

def seq_to_codes(seq):
    """ACGT string -> uint8 codes 0-3; anything else is 4."""
    return _BASE_LUT[np.frombuffer(seq.encode("ascii"), dtype=np.uint8)]


def revcomp_codes(codes):
    out = codes[::-1].copy()
    ok = out < 4
    out[ok] = 3 - out[ok]
    return out


def kmer_codes(codes, k):
    """Base-4 code of every k-mer in a code array (len n-k+1); -1 if it has a non-ACGT."""
    if codes.size < k:
        return np.empty(0, dtype=np.int64)
    win = np.lib.stride_tricks.sliding_window_view(codes, k)
    weights = 4 ** np.arange(k - 1, -1, -1, dtype=np.int64)
    out = win.astype(np.int64) @ weights
    out[(win > 3).any(axis=1)] = -1
    return out


def code_to_kmer(code, k):
    return "".join("ACGT"[(int(code) >> (2 * (k - i - 1))) & 3] for i in range(k))


def kmer_to_code(kmer):
    code = 0
    for b in kmer:
        code = (code << 2) | "ACGT".index(b)
    return code


class RefCache:
    """One contig resident at a time: its sequence and its base codes."""

    def __init__(self, fasta_path):
        import pysam
        self.fasta = pysam.FastaFile(fasta_path)
        self.name = None
        self.seq = ""
        self.codes = np.empty(0, dtype=np.uint8)

    def get(self, name):
        if name != self.name:
            try:
                self.seq = self.fasta.fetch(name).upper()
            except (ValueError, KeyError):
                sys.exit(f"Error: contig {name} is in the BAM but not in the reference "
                         f"FASTA. The k-mer context cannot be determined without it.")
            self.codes = seq_to_codes(self.seq)
            self.name = name
        return self.seq, self.codes

    def close(self):
        self.fasta.close()


# ---------------------------------------------------------------------------
# Which reads, and which observations, are usable
# ---------------------------------------------------------------------------

def read_half(read_id):
    """Deterministic split of reads into halves 'A' and 'B' by read name.

    The control sample is scored against a model built from itself. Doing that
    with the same reads is circular - every read pulls the model toward its own
    level - so the model takes half A and the control is scored on half B. The
    split depends only on the read name, so both backends split identically.
    """
    return "A" if zlib.crc32(read_id.encode()) & 1 == 0 else "B"


def ref_error_mask(sam, ref_codes):
    """Bool array over [reference_start, reference_end): True where this read does
    not match the reference - mismatch, deleted base, or a base flanking an insertion.

    Read from CIGAR + query sequence + reference, so no MD tag is needed (the
    toolkit's minimap2 call does not write one). Returns None if the record carries
    no sequence.
    """
    seq = sam.query_sequence
    if not seq or sam.cigartuples is None:
        return None
    q = seq_to_codes(seq)
    start = sam.reference_start
    err = np.zeros(sam.reference_end - start, dtype=bool)
    r = 0
    qi = 0
    for op, ln in sam.cigartuples:
        if op in (0, 7, 8):                       # M, =, X
            err[r:r + ln] |= q[qi:qi + ln] != ref_codes[start + r:start + r + ln]
            r += ln
            qi += ln
        elif op in (2, 3):                        # D, N
            err[r:r + ln] = True
            r += ln
        elif op == 1:                             # I
            if r > 0:
                err[r - 1] = True
            if r < err.size:
                err[r] = True
            qi += ln
        elif op == 4:                             # S
            qi += ln
    return err


def clean_positions(err, start, pos, flank):
    """True for positions with no alignment error within +/- flank, entirely inside
    the aligned span. flank = k-1 covers every k-mer that could contain the position
    whatever base of the k-mer the level is registered to, so this does not depend
    on either backend's registration.
    """
    if err is None:
        return np.zeros(pos.size, dtype=bool)
    cs = np.concatenate(([0], np.cumsum(err)))
    lo = pos - start - flank
    hi = pos - start + flank + 1
    inside = (lo >= 0) & (hi <= err.size)
    out = np.zeros(pos.size, dtype=bool)
    out[inside] = (cs[hi[inside]] - cs[lo[inside]]) == 0
    return out


class ReadLevels:
    """One read's observations. pos is 0-based reference coordinates; kmer is -1
    where the context holds a non-ACGT base; clean is all False unless requested."""
    __slots__ = ("ref", "start", "read_id", "pos", "level", "kmer", "clean")

    def __init__(self, ref, start, read_id, pos, level, kmer, clean):
        self.ref = ref
        self.start = start
        self.read_id = read_id
        self.pos = pos
        self.level = level
        self.kmer = kmer
        self.clean = clean


# ---------------------------------------------------------------------------
# Backends
# ---------------------------------------------------------------------------

class _Source:
    """Common filtering and bookkeeping for the two backends."""

    def __init__(self, bam, genome, strand, half, need_clean):
        self.bam = bam
        self.genome = genome
        self.want_reverse = None if strand is None else (strand == "rev")
        self.half = half
        self.need_clean = need_clean
        self.refs = RefCache(genome)
        self.counts = {"candidate reads": 0, "reads used": 0, "reads failed": 0,
                       "positions without a level": 0}

    def wanted(self, sam):
        if self.want_reverse is not None and bool(sam.is_reverse) != self.want_reverse:
            return False
        if self.half != "all" and read_half(sam.query_name) != self.half:
            return False
        return True

    def clean_for(self, sam, pos):
        if not self.need_clean:
            return np.zeros(pos.size, dtype=bool)
        _, codes = self.refs.get(sam.reference_name)
        err = ref_error_mask(sam, codes)
        return clean_positions(err, sam.reference_start, pos, self.kmer_len - 1)

    def report(self):
        for key, n in self.counts.items():
            sys.stderr.write(f"{key}: {n}\n")
        # A scoring or model run that silently used nothing looks like a sample with
        # no signal. Reads went in, so something must come out.
        if self.counts["candidate reads"] > 0 and self.counts["reads used"] == 0:
            sys.exit(f"Error: {self.counts['candidate reads']} reads were selected from "
                     f"{self.bam} but none yielded signal levels. See the counts above.")


class Uncalled4Levels(_Source):
    """dtw.current / seq.kmer / seq.pos from an Uncalled4 BAM, via uncalled4's own
    decoder - for the reasons spelled out in perbase_signal_deviation.py: the k-mer
    shift and the reverse-strand coordinate flip are its code, not ours.

    The level is dtw.current: mean normalized current of the samples DTW assigned to
    the reference position (aln.hpp), in the pore model's normalized units.
    """
    name = "uncalled4"

    def __init__(self, bam, genome, strand=None, half="all", need_clean=False, **_):
        super().__init__(bam, genome, strand, half, need_clean)
        try:
            from uncalled4 import Config, Tracks
        except ImportError as exc:
            sys.exit(f"Error: could not import uncalled4.\n"
                     f"  {type(exc).__name__}: {exc}\n  python:  {sys.executable}")
        conf = Config()
        conf.tracks.io.bam_in = [bam]
        conf.tracks.ref = genome
        conf.tracks.layers = ["dtw"]
        conf.read_index.load_signal = False
        self.tracks = Tracks(conf=conf)
        self.bam_in = self.tracks.bam_in
        if self.bam_in is None:
            sys.exit(f"Error: uncalled4 did not open {bam} as a BAM input.")
        model = self.tracks.model
        if model is None:
            sys.exit("Error: uncalled4 loaded no pore model for this BAM; the k-mer "
                     "length is unknown.")
        self.kmer_len = int(model.K)
        shift = int(model.PRMS.shift)
        self.context = (shift, self.kmer_len - shift - 1)
        # The model file names k-mers by decoding uncalled4's integers with OUR
        # base-4 layout. Scoring is unaffected if the layouts differ - the same
        # decode/encode round trip is applied on both sides - but the k-mer strings
        # in the model file would then be wrong, so say so. (A list, not a str:
        # str_to_kmer iterates a bare string character by character.)
        probe = ("ACGT" * 8)[:self.kmer_len - 1] + "T"
        try:
            got = int(np.asarray(model.str_to_kmer([probe])).ravel()[0])
        except Exception as exc:
            got = f"{type(exc).__name__}: {exc}"
        if got != kmer_to_code(probe):
            sys.stderr.write(
                f"Warning: uncalled4 encodes {probe} as {got}, expected "
                f"{kmer_to_code(probe)}. K-mer names in the model file may not be the "
                f"true sequences; calls are unaffected.\n")

    @staticmethod
    def _np(x, dtype):
        return np.asarray(x.to_numpy() if hasattr(x, "to_numpy") else x, dtype=dtype)

    def __iter__(self):
        for sam in self.bam_in.iter_sam():
            if not self.wanted(sam):
                continue
            self.counts["candidate reads"] += 1
            aln = self.bam_in.sam_to_aln(sam, load_moves=False)
            if aln is None:
                self.counts["reads failed"] += 1
                continue
            pos = self._np(aln.seq.pos, np.int64)
            level = self._np(aln.dtw.current, np.float64)
            kmer = self._np(aln.seq.kmer, np.int64)
            if not (pos.size == level.size == kmer.size):
                sys.exit(f"Error: read {sam.query_name} has {pos.size} positions, "
                         f"{level.size} levels and {kmer.size} k-mers. Refusing to "
                         f"guess how they line up.")
            keep = ~np.isnan(level)
            self.counts["positions without a level"] += int((~keep).sum())
            pos, level, kmer = pos[keep], level[keep], kmer[keep]
            self.counts["reads used"] += 1
            yield ReadLevels(sam.reference_name, sam.reference_start, sam.query_name,
                             pos, level, kmer, self.clean_for(sam, pos))

    def close(self):
        self.tracks.close()
        self.refs.close()


class RemoraLevels(_Source):
    """Per-base trimmed-mean level from pod5 + BAM through the Remora API.

    io.Read.from_pod5_and_alignment -> (optional) set_refine_signal_mapping ->
    compute_per_base_metric("dwell_trimmean_trimsd"), reference-anchored. The level
    is Remora's normalized signal, trimmed by one sample at each end of the base.

    The BAM needs the dorado move table (mv, ts) and must be coordinate-sorted and
    indexed. The Uncalled4 BAM qualifies - it keeps those tags - and using it means
    both backends score exactly the same reads.

    Without --levels there is no signal-mapping refinement: boundaries are the
    basecaller's moves at stride resolution, and the base a level belongs to is not
    known. With a level table (ONT kmer_models, e.g.
    dna_r10.4.1_e8.2_400bps/9mer_levels_v1.txt) the mapping is rescaled and refined.

    Either way the k-mer context is a CENTRED 9-mer (4 bases each side) unless
    --kmer-context says otherwise: the R10 pore has two reader heads, so the current
    at a base depends on sequence on both sides of it rather than on one dominant
    position, and the level table's own dominant position is not used for the context.
    """
    name = "remora"
    BATCH = 2000

    def __init__(self, bam, genome, strand=None, half="all", need_clean=False,
                 pod5=None, levels=None, context=None, **_):
        super().__init__(bam, genome, strand, half, need_clean)
        if not pod5:
            sys.exit("Error: the remora backend needs --pod5 (raw signal).")
        try:
            import pod5 as pod5_lib
            from remora import io, refine_signal_map, util, RemoraError
        except ImportError as exc:
            sys.exit(f"Error: could not import remora/pod5.\n"
                     f"  {type(exc).__name__}: {exc}\n  python:  {sys.executable}\n"
                     f"Install with: pip install ont-remora pod5")
        self.io = io
        self.revcomp = util.revcomp
        self.RemoraError = RemoraError
        self.pod5_dr = pod5_lib.DatasetReader(pod5, recursive=True)
        if levels:
            self.refiner = refine_signal_map.SigMapRefiner(
                kmer_model_filename=levels, do_rough_rescale=True, scale_iters=0,
                do_fix_guage=True)
        else:
            sys.stderr.write("Warning: no k-mer level table (--levels); signal mapping "
                             "is NOT refined and levels use the basecaller's moves.\n")
            self.refiner = None
        self.context = tuple(context) if context else (4, 4)
        self.kmer_len = sum(self.context) + 1
        self.counts["reads missing from pod5"] = 0

    def _levels(self, sam, p5):
        io_read = self.io.Read.from_pod5_and_alignment(p5, sam)
        seq, codes = self.refs.get(sam.reference_name)
        start, end = sam.reference_start, sam.reference_end
        if io_read.ref_seq is None:
            # No MD tag: take the reference span from the FASTA instead.
            span = seq[start:end]
            io_read.ref_seq = self.revcomp(span) if sam.is_reverse else span
            io_read.compute_ref_to_signal()
        if self.refiner is not None:
            io_read.set_refine_signal_mapping(self.refiner, ref_mapping=True)
        level = np.asarray(
            io_read.compute_per_base_metric("dwell_trimmean_trimsd")["trimmean"],
            dtype=np.float64)
        n = end - start
        if level.size != n:
            raise self.RemoraError("per-base metric length != reference span")

        # Metrics and k-mers are in read orientation; index i is the i-th base the
        # pore saw. Pad the reference span by the context, clipping at contig ends.
        before, after = self.context
        left, right = (after, before) if sam.is_reverse else (before, after)
        lo, hi = start - left, end + right
        ext = np.full(hi - lo, 4, dtype=np.uint8)
        a, b = max(lo, 0), min(hi, codes.size)
        ext[a - lo:b - lo] = codes[a:b]
        if sam.is_reverse:
            ext = revcomp_codes(ext)
            pos = np.arange(end - 1, start - 1, -1, dtype=np.int64)
        else:
            pos = np.arange(start, end, dtype=np.int64)
        return pos, level, kmer_codes(ext, self.kmer_len)

    def _batch(self, sams):
        ids = {}
        for sam in sams:
            ids[sam.get_tag("pi") if sam.has_tag("pi") else sam.query_name] = None
        for p5 in self.pod5_dr.reads(selection=list(ids), missing_ok=True,
                                     preload=["samples"]):
            ids[str(p5.read_id)] = p5
        for sam in sams:
            p5 = ids[sam.get_tag("pi") if sam.has_tag("pi") else sam.query_name]
            if p5 is None:
                self.counts["reads missing from pod5"] += 1
                continue
            try:
                pos, level, kmer = self._levels(sam, p5)
            except self.RemoraError:
                self.counts["reads failed"] += 1
                continue
            keep = np.isfinite(level)
            self.counts["positions without a level"] += int((~keep).sum())
            pos, level, kmer = pos[keep], level[keep], kmer[keep]
            self.counts["reads used"] += 1
            yield ReadLevels(sam.reference_name, sam.reference_start, sam.query_name,
                             pos, level, kmer, self.clean_for(sam, pos))

    def __iter__(self):
        import pysam
        batch = []
        with pysam.AlignmentFile(self.bam) as bam_fh:
            for sam in bam_fh.fetch():
                if sam.is_unmapped or sam.is_secondary or sam.is_supplementary:
                    continue
                if not self.wanted(sam):
                    continue
                self.counts["candidate reads"] += 1
                batch.append(sam)
                if len(batch) >= self.BATCH:
                    yield from self._batch(batch)
                    batch = []
            yield from self._batch(batch)

    def close(self):
        self.pod5_dr.close()
        self.refs.close()


def add_backend_args(p, multi=False):
    """multi: accept several BAMs (and, for remora, as many pod5 paths) to pool."""
    nargs = "+" if multi else None
    p.add_argument("--backend", required=True, choices=BACKENDS,
                   help="Where signal levels come from: the Uncalled4 BAM's DTW tags, "
                        "or pod5 + move table through the Remora API.")
    p.add_argument("-i", "--bam", required=True, nargs=nargs,
                   help="Coordinate-sorted, indexed BAM. uncalled4: an Uncalled4 BAM. "
                        "remora: any aligned BAM that kept the dorado mv/ts tags "
                        "(the Uncalled4 BAM does).")
    p.add_argument("-g", "--genome", required=True, help="Reference genome FASTA")
    p.add_argument("--pod5", nargs=nargs,
                   help="remora: pod5 file or directory (searched recursively)"
                        + ("; one per BAM, in the same order" if multi else ""))
    p.add_argument("--levels", help="remora: k-mer level table for signal mapping "
                                    "refinement (ONT kmer_models). Omit to use raw moves.")
    p.add_argument("--kmer-context", type=int, nargs=2, metavar=("BEFORE", "AFTER"),
                   help="remora: bases before/after the level's base that define the "
                        "k-mer (default: 4 4, a centred 9-mer)")
    p.add_argument("--half", default="all", choices=["all", "A", "B"],
                   help="Use all reads, or one half of a fixed split by read name "
                        "(model on A, score the same sample on B)")


def open_source(args, strand, need_clean, bam=None, pod5=None):
    """bam/pod5 override args.bam/args.pod5, for scripts that take several."""
    cls = Uncalled4Levels if args.backend == "uncalled4" else RemoraLevels
    return cls(bam or args.bam, args.genome, strand=strand, half=args.half,
               need_clean=need_clean, pod5=pod5 or args.pod5, levels=args.levels,
               context=args.kmer_context)


# ---------------------------------------------------------------------------
# Null distribution of the calling statistic, for a false-positive-rate threshold
# ---------------------------------------------------------------------------

HIST_WIDTH = 0.01
HIST_BINS = 10000           # statistic 0 .. 100, plus one overflow bin


def hist_add(hist, stat):
    np.add.at(hist, np.minimum((stat / HIST_WIDTH).astype(np.int64), HIST_BINS), 1)


def write_hist(path, hist, meta):
    with open(path, "w") as out:
        for key, value in meta.items():
            out.write(f"#{key}={value}\n")
        out.write("bin_lower\tcount\n")
        for i in np.flatnonzero(hist):
            out.write(f"{i * HIST_WIDTH:.2f}\t{int(hist[i])}\n")


def cutoff_from_hists(paths, fpr, expect_meta):
    """Smallest statistic cutoff that at most a fraction `fpr` of the null
    observations exceed. The null is the held-out control scored against the model,
    pooled over every histogram given (both strands; every member of a pooled model).
    """
    hist = np.zeros(HIST_BINS + 1, dtype=np.int64)
    for path in paths:
        with open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    key, _, value = line[1:].rstrip("\n").partition("=")
                    if key in expect_meta and str(expect_meta[key]) != value:
                        sys.exit(f"Error: null histogram {path} has {key}={value}, "
                                 f"this run has {key}={expect_meta[key]}.")
                elif not line.startswith("bin_lower"):
                    lower, count = line.split("\t")
                    hist[int(round(float(lower) / HIST_WIDTH))] += int(count)
    total = int(hist.sum())
    if total == 0:
        sys.exit("Error: the null histograms hold no observations; a false-positive "
                 "rate cannot be estimated from an empty held-out control.")
    tail = np.cumsum(hist[::-1])[::-1]          # observations in bin i or above
    first = int(np.argmax(tail <= fpr * total)) if (tail <= fpr * total).any() \
        else HIST_BINS + 1
    return first * HIST_WIDTH, total


# ---------------------------------------------------------------------------
# K-mer model
# ---------------------------------------------------------------------------

def _group_bounds(sorted_keys):
    starts = np.concatenate(([0], np.flatnonzero(np.diff(sorted_keys)) + 1))
    counts = np.diff(np.concatenate((starts, [sorted_keys.size])))
    return starts, counts


def _group_median(sorted_vals, starts, counts):
    lo = starts + (counts - 1) // 2
    hi = starts + counts // 2
    return 0.5 * (sorted_vals[lo] + sorted_vals[hi])


class KmerModelBuilder:
    """Accumulates (k-mer, level) observations and reduces them to n / mean / SD.

    max_obs caps what is held per k-mer with an exact uniform sample: every
    observation gets a random priority and the max_obs smallest per k-mer survive
    each compaction. That bounds memory on deep or genome-scale data without
    favouring the reads that happen to come first in a coordinate-sorted BAM.
    """

    COMPACT_AT = 20_000_000

    def __init__(self, max_obs=0, seed=1):
        self.max_obs = max_obs
        self.rng = np.random.default_rng(seed)
        self.parts = []
        self.buffered = 0
        self.seen = 0

    def add(self, kmer, level):
        if kmer.size == 0:
            return
        prio = self.rng.random(kmer.size) if self.max_obs else None
        self.parts.append((kmer, level, prio))
        self.buffered += kmer.size
        self.seen += kmer.size
        if self.max_obs and self.buffered >= self.COMPACT_AT:
            self._compact()

    def _concat(self):
        kmer = np.concatenate([p[0] for p in self.parts])
        level = np.concatenate([p[1] for p in self.parts])
        prio = np.concatenate([p[2] for p in self.parts]) if self.max_obs else None
        return kmer, level, prio

    def _compact(self):
        kmer, level, prio = self._concat()
        order = np.lexsort((prio, kmer))
        kmer, level, prio = kmer[order], level[order], prio[order]
        starts, counts = _group_bounds(kmer)
        rank = np.arange(kmer.size) - np.repeat(starts, counts)
        keep = rank < self.max_obs
        self.parts = [(kmer[keep], level[keep], prio[keep])]
        self.buffered = int(keep.sum())

    def finalize(self, mad_k):
        """-> dict of arrays: kmer, n, n_outlier, mean, sd, median, mad (sorted by kmer)."""
        if not self.parts:
            return {key: np.empty(0) for key in
                    ("kmer", "n", "n_outlier", "mean", "sd", "median", "mad")}
        if self.max_obs:
            self._compact()
        kmer, level, _ = self._concat()
        order = np.lexsort((level, kmer))
        kmer, level = kmer[order], level[order]
        starts, counts = _group_bounds(kmer)
        median = _group_median(level, starts, counts)

        absdev = np.abs(level - np.repeat(median, counts))
        mad = _group_median(absdev[np.lexsort((absdev, kmer))], starts, counts)
        keep = absdev <= mad_k * np.repeat(mad, counts)

        n_all = counts
        kmer_k, level_k = kmer[keep], level[keep]
        starts_k, n = _group_bounds(kmer_k)        # every k-mer keeps its median
        mean = np.add.reduceat(level_k, starts_k) / n
        ss = np.add.reduceat((level_k - np.repeat(mean, n)) ** 2, starts_k)
        with np.errstate(invalid="ignore", divide="ignore"):
            sd = np.sqrt(ss / (n - 1))
        return {"kmer": kmer_k[starts_k], "n": n, "n_outlier": n_all - n,
                "mean": mean, "sd": sd, "median": median, "mad": mad}


MODEL_COLUMNS = ("kmer", "n", "n_outlier", "mean", "sd", "median", "mad")


def write_model(path, meta, table):
    k = int(meta["kmer_len"])
    with open(path, "w") as out:
        for key, value in meta.items():
            out.write(f"#{key}={value}\n")
        out.write("\t".join(MODEL_COLUMNS) + "\n")
        for i in range(table["kmer"].size):
            out.write(f"{code_to_kmer(table['kmer'][i], k)}\t{int(table['n'][i])}"
                      f"\t{int(table['n_outlier'][i])}\t{table['mean'][i]:.6f}"
                      f"\t{table['sd'][i]:.6f}\t{table['median'][i]:.6f}"
                      f"\t{table['mad'][i]:.6f}\n")


def load_model(path):
    """-> (meta dict, kmer codes sorted, n, mean, sd)."""
    meta = {}
    kmer, n, mean, sd = [], [], [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                key, _, value = line[1:].rstrip("\n").partition("=")
                meta[key] = value
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] == "kmer":
                continue
            kmer.append(kmer_to_code(f[0]))
            n.append(int(f[1]))
            mean.append(float(f[3]))
            sd.append(float(f[4]))
    kmer = np.array(kmer, dtype=np.int64)
    order = np.argsort(kmer)
    return (meta, kmer[order], np.array(n, dtype=np.int64)[order],
            np.array(mean)[order], np.array(sd)[order])


# ---------------------------------------------------------------------------
# The test
# ---------------------------------------------------------------------------

def z_scores(kmer, level, m_kmer, m_n, m_mean, m_sd, min_obs):
    """Signed z against the model; NaN where the k-mer is not modelled well enough."""
    z = np.full(level.size, np.nan)
    if m_kmer.size == 0:
        return z
    idx = np.minimum(np.searchsorted(m_kmer, kmer), m_kmer.size - 1)
    ok = (m_kmer[idx] == kmer) & (m_n[idx] >= min_obs) & (m_sd[idx] > 0)
    z[ok] = (level[ok] - m_mean[idx[ok]]) / m_sd[idx[ok]]
    return z


try:
    from scipy.special import erfc as _erfc
except ImportError:
    _erfc_py = np.frompyfunc(math.erfc, 1, 1)

    def _erfc(x):
        return _erfc_py(x).astype(np.float64)


def two_sided_p(z):
    """P(|Z| >= |z|) under a standard normal - Tombo's norm.cdf(-|z|) * 2."""
    return _erfc(np.abs(z) / math.sqrt(2.0))


def z_threshold(pval):
    """|z| at which the two-sided p-value equals pval."""
    return NormalDist().inv_cdf(1.0 - pval / 2.0)


def fisher_window(pos, p, lag):
    """Fisher's method over reference positions pos-lag..pos+lag of ONE read.

    Tombo's calc_window_fishers_method, except that a window missing some positions
    (read end, deletion, unmodelled k-mer) combines the ones present with matching
    degrees of freedom instead of returning NaN - Tombo's own
    calc_vectorized_fm_pvals(filter_nan=True) behaviour. chi2 with 2m degrees of
    freedom has the closed-form survival function exp(-s) * sum_{i<m} s^i / i!
    for s = -sum(log p), so no scipy is needed.
    """
    if pos.size == 0:
        return p
    lo = int(pos.min())
    span = int(pos.max()) - lo + 1
    logp = np.zeros(span)
    present = np.zeros(span)
    logp[pos - lo] = np.log(np.maximum(p, SMALLEST_P))
    present[pos - lo] = 1.0
    kernel = np.ones(2 * lag + 1)
    s = -np.convolve(logp, kernel, mode="same")[pos - lo]
    m = np.convolve(present, kernel, mode="same")[pos - lo].round().astype(np.int64)
    term = np.ones(pos.size)
    total = np.zeros(pos.size)
    for i in range(2 * lag + 1):
        total += np.where(i < m, term, 0.0)
        term = term * s / (i + 1)
    return np.minimum(np.exp(-s) * total, 1.0)


# ---------------------------------------------------------------------------
# Per-position accumulation
# ---------------------------------------------------------------------------

class PositionAccumulator:
    """Streams per-read calls into per-position totals for a coordinate-sorted BAM.

    Dense blocks of BLOCK positions are finalised once a read starts past them, the
    same watermark logic as perbase_signal_deviation.py's WindowAccumulator, and for
    the same reason an observation below the watermark is an error, not a drop.
    Memory is 32 bytes per position in the open blocks, independent of coverage.
    """
    BLOCK = 1 << 16

    def __init__(self, emit):
        self.emit = emit            # emit(ref, positions, cov, n_mod, sum_z, sum_abs_z)
        self.ref = None
        self.blocks = {}
        self.watermark = 0

    def _flush(self, below=None):
        for b in sorted(k for k in self.blocks if below is None or k < below):
            a = self.blocks.pop(b)
            idx = np.flatnonzero(a[0])
            self.emit(self.ref, idx + b * self.BLOCK, a[0][idx], a[1][idx],
                      a[2][idx], a[3][idx])
        if below is not None:
            self.watermark = max(self.watermark, below)

    def set_ref(self, ref):
        if ref != self.ref:
            self._flush()
            self.ref = ref
            self.watermark = 0

    def add(self, pos, z, called, read_start):
        if pos.size:
            blk = pos // self.BLOCK
            if int(blk.min()) < self.watermark:
                raise RuntimeError(
                    f"Observation at {self.ref}:{int(pos.min())} is below the "
                    f"finalised watermark. The BAM is not coordinate-sorted, or a "
                    f"read reports positions before its own alignment start.")
            for b in np.unique(blk):
                m = blk == b
                a = self.blocks.get(int(b))
                if a is None:
                    a = np.zeros((4, self.BLOCK))
                    self.blocks[int(b)] = a
                off = pos[m] - int(b) * self.BLOCK
                np.add.at(a[0], off, 1.0)
                np.add.at(a[1], off, called[m].astype(np.float64))
                np.add.at(a[2], off, z[m])
                np.add.at(a[3], off, np.abs(z[m]))
        self._flush(below=read_start // self.BLOCK)

    def close(self):
        self._flush()


def open_out(path):
    return gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w")
