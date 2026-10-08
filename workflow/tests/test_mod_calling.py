#!/usr/bin/env python3
"""
test_mod_calling.py - Tests for direct modification calling (phases 2c/3c) that need
neither uncalled4 nor remora nor any real data.

Everything downstream of "a read's (position, level, k-mer)" is shared by the two
backends, so it is tested here on simulated reads fed through a stand-in source:
the k-mer model, the MAD filter, the error mask, the z/p test, Fisher's method, the
per-position accumulator, both command-line scripts end to end, and mod_reactivity.

What this does NOT test is the two backends themselves - that dtw.current / seq.kmer
decode as expected from a real Uncalled4 BAM, and that Remora's per-base metric
lines up with reference positions. Those need real data; see the README.

Usage:
    python3 workflow/tests/test_mod_calling.py
"""

import os
import sys
import gzip
import math
import tempfile
import types

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "scripts"))
import mod_calling_common as mc
import kmer_signal_model
import perbase_mod_calls
import mod_reactivity

K = 9
SHIFT = 6
RNG = np.random.default_rng(7)
GENOME = {"oligoA": "".join(RNG.choice(list("ACGT"), 120)),
          "oligoB": "".join(RNG.choice(list("ACGT"), 90))}
LEVEL = RNG.normal(0, 1, 4 ** K)          # true level of every k-mer
NOISE = 0.1


class FakeSam:
    def __init__(self, query_sequence, cigartuples, reference_start, reference_end):
        self.query_sequence = query_sequence
        self.cigartuples = cigartuples
        self.reference_start = reference_start
        self.reference_end = reference_end


def simulate(ref, n_reads, reverse, modified=None, prefix="r", rng=None):
    """Full-length reads of one strand. `modified` maps 0-based position -> (rate,
    shift): that fraction of reads has its level moved by `shift` there."""
    rng = rng or np.random.default_rng(1)
    codes = mc.seq_to_codes(GENOME[ref])
    before, after = SHIFT, K - SHIFT - 1
    reads = []
    for i in range(n_reads):
        if reverse:
            # read orientation: k-mer j of the reverse complement is registered to
            # its base j+before, which is reference position L-1-(j+before)
            km = mc.kmer_codes(mc.revcomp_codes(codes), K)
            pos = codes.size - 1 - (np.arange(km.size) + before)
        else:
            km = mc.kmer_codes(codes, K)
            pos = np.arange(before, codes.size - after)
        level = LEVEL[km] + rng.normal(0, NOISE, km.size)
        for p, (rate, shift) in (modified or {}).items():
            if rng.random() < rate:
                level[pos == p] += shift
        reads.append(mc.ReadLevels(ref, 0, f"{prefix}{i}", pos.astype(np.int64),
                                   level, km.astype(np.int64),
                                   np.ones(km.size, dtype=bool)))
    return reads


class FakeSource:
    """Stands in for Uncalled4Levels / RemoraLevels."""
    kmer_len = K
    context = (SHIFT, K - SHIFT - 1)

    def __init__(self, reads, fasta, half):
        self.reads = reads
        self.half = half
        self.refs = mc.RefCache(fasta)
        self.counts = {"candidate reads": 0}

    def __iter__(self):
        for r in self.reads:
            if self.half == "all" or mc.read_half(r.read_id) == self.half:
                self.counts["candidate reads"] += 1
                yield r

    def report(self):
        pass

    def close(self):
        self.refs.close()


def run(module, argv, reads_for):
    """Run a script's main() with open_source replaced by simulated reads.
    reads_for(strand, bam) -> reads; bam is set when the script pools several."""
    def fake_open(args, strand, need_clean, bam=None, pod5=None):
        return FakeSource(reads_for(strand, bam), args.genome, args.half)
    old_open, old_argv = module.open_source, sys.argv
    module.open_source, sys.argv = fake_open, [module.__name__] + argv
    try:
        module.main()
    finally:
        module.open_source, sys.argv = old_open, old_argv


def read_table(path):
    with gzip.open(path, "rt") as fh:
        return [line.rstrip("\n").split("\t") for line in fh]


def test_encoding():
    assert mc.kmer_to_code("AAAAAAAAT") == 3
    assert mc.code_to_kmer(mc.kmer_to_code("ACGTTGCAA"), 9) == "ACGTTGCAA"
    codes = mc.seq_to_codes("ACGTNACGT")
    assert list(mc.kmer_codes(codes, 4)) == [27, -1, -1, -1, -1, 27]
    assert list(mc.revcomp_codes(mc.seq_to_codes("AACGN"))) == [4, 1, 2, 3, 3]


def test_error_mask():
    ref = mc.seq_to_codes("ACGTACGTACGTACGTACGT")
    # 5M 1I 4M 2D 5M with a mismatch at reference position 3
    sam = FakeSam("ACGAA" + "T" + "CGTA" + "TACGT",
                  [(0, 5), (1, 1), (0, 4), (2, 2), (0, 5)], 0, 16)
    err = mc.ref_error_mask(sam, ref)
    assert list(np.flatnonzero(err)) == [3, 4, 5, 9, 10], np.flatnonzero(err)
    clean = mc.clean_positions(err, 0, np.arange(16), 1)
    # within +/-1 of an error, or within 1 of either alignment end -> not clean
    assert list(np.flatnonzero(clean)) == [1, 7, 12, 13, 14], np.flatnonzero(clean)
    assert mc.ref_error_mask(FakeSam(None, [(0, 5)], 0, 5), ref) is None


def test_model_and_mad():
    b = mc.KmerModelBuilder(max_obs=0)
    vals = RNG.normal(5.0, 0.2, 2000)
    b.add(np.full(2000, 11, dtype=np.int64), vals)
    b.add(np.array([11, 11, 42], dtype=np.int64), np.array([500.0, -500.0, 1.0]))
    t = b.finalize(15.0)
    assert list(t["kmer"]) == [11, 42]
    assert list(t["n"]) == [2000, 1] and list(t["n_outlier"]) == [2, 0]
    assert abs(t["mean"][0] - vals.mean()) < 1e-9
    assert abs(t["sd"][0] - vals.std(ddof=1)) < 1e-9
    assert np.isnan(t["sd"][1])            # one observation: no SD, never scored

    # the per-k-mer cap is a uniform sample, not "the first N"
    b = mc.KmerModelBuilder(max_obs=500, seed=3)
    for chunk in np.array_split(np.arange(100000, dtype=np.float64), 50):
        b.add(np.zeros(chunk.size, dtype=np.int64), chunk)
    t = b.finalize(1e9)
    assert t["n"][0] == 500
    assert 40000 < t["mean"][0] < 60000, t["mean"][0]


def test_statistics():
    z = np.array([0.0, 1.0, -2.3263478740408408, 5.0])
    p = mc.two_sided_p(z)
    assert abs(p[0] - 1.0) < 1e-12 and abs(p[2] - 0.02) < 1e-9
    assert abs(mc.z_threshold(0.02) - 2.3263478740408408) < 1e-9

    # Fisher: chi2 survival with 2m dof, checked against the closed form by hand
    pos = np.array([10, 11, 12, 20])
    pv = np.array([0.01, 0.5, 0.2, 0.03])
    out = mc.fisher_window(pos, pv, 1)

    def chi2_sf(ps):
        s = -sum(math.log(x) for x in ps)
        return math.exp(-s) * sum(s ** i / math.factorial(i) for i in range(len(ps)))
    expect = [chi2_sf([0.01, 0.5]), chi2_sf([0.01, 0.5, 0.2]), chi2_sf([0.5, 0.2]), 0.03]
    assert np.allclose(out, expect), (out, expect)


def test_accumulator():
    got = {}

    def emit(ref, positions, cov, n_mod, sum_z, sum_abs):
        for i in range(positions.size):
            got[(ref, int(positions[i]))] = (cov[i], n_mod[i], sum_z[i], sum_abs[i])
    acc = mc.PositionAccumulator(emit)
    B = acc.BLOCK
    acc.set_ref("c1")
    acc.add(np.array([B - 1, B]), np.array([1.0, -3.0]), np.array([False, True]), B - 1)
    acc.add(np.array([B, B + 5]), np.array([-1.0, 2.0]), np.array([True, False]), B)
    assert (("c1", B - 1) in got) and (("c1", B) not in got)     # block 0 finalised
    acc.set_ref("c2")
    acc.add(np.array([3]), np.array([4.0]), np.array([True]), 3)
    acc.close()
    assert got[("c1", B)] == (2.0, 2.0, -4.0, 4.0)
    assert got[("c2", 3)] == (1.0, 1.0, 4.0, 4.0)
    try:
        acc.set_ref("c3")
        acc.add(np.array([2 * B]), np.array([0.0]), np.array([False]), 2 * B)
        acc.add(np.array([5]), np.array([0.0]), np.array([False]), 2 * B)
    except RuntimeError:
        pass
    else:
        raise AssertionError("an observation below the watermark must be an error")


def test_end_to_end():
    tmp = tempfile.mkdtemp(prefix="mod_calling_test_")
    fasta = os.path.join(tmp, "genome.fa")
    with open(fasta, "w") as fh:
        for name, seq in GENOME.items():
            fh.write(f">{name}\n{seq}\n")
    import pysam
    pysam.faidx(fasta)
    bam = os.path.join(tmp, "placeholder.bam")
    open(bam, "w").close()
    common = ["--backend", "uncalled4", "-i", bam, "-g", fasta]

    # Treated: position 40 modified in 60% of forward reads, 70 in 25%.
    mods = {40: (0.60, 1.0), 70: (0.25, -1.0)}
    control = {s: simulate("oligoA", 4000, s == "rev", prefix=f"c{s}",
                           rng=np.random.default_rng(11 + (s == "rev")))
               + simulate("oligoB", 4000, s == "rev", prefix=f"cb{s}",
                          rng=np.random.default_rng(13 + (s == "rev")))
               for s in ("for", "rev")}
    treated = {s: simulate("oligoA", 3000, s == "rev", prefix=f"t{s}",
                           modified=mods if s == "for" else None,
                           rng=np.random.default_rng(21 + (s == "rev")))
               + simulate("oligoB", 3000, s == "rev", prefix=f"tb{s}",
                          rng=np.random.default_rng(23 + (s == "rev")))
               for s in ("for", "rev")}

    model = os.path.join(tmp, "model.tsv")
    run(kmer_signal_model, common + ["--half", "A", "-o", model],
        lambda strand, bam: control["for"] + control["rev"])
    meta, m_kmer, m_n, m_mean, m_sd = mc.load_model(model)
    assert meta["backend"] == "uncalled4" and meta["kmer_len"] == "9"
    n_expect = 2 * sum(len(s) - K + 1 for s in GENOME.values())
    assert m_kmer.size == n_expect, (m_kmer.size, n_expect)
    assert np.allclose(m_mean, LEVEL[m_kmer], atol=0.02)
    assert np.allclose(m_sd, NOISE, rtol=0.15)

    out = {}
    for name, reads, half in (("ctrl", control, "B"), ("trt", treated, "all")):
        for strand in ("for", "rev"):
            path = os.path.join(tmp, f"{name}_{strand}.txt.gz")
            per_read = os.path.join(tmp, f"{name}_{strand}_reads.tsv.gz")
            run(perbase_mod_calls,
                common + ["-m", model, "-s", strand, "--half", half, "-o", path,
                          "--per-read", per_read],
                lambda s, bam, reads=reads: reads[s])
            out[(name, strand)] = read_table(path)

    # Held-out control: the false-positive rate is the p-value threshold.
    ctrl = out[("ctrl", "for")]
    frac = np.array([float(r[4]) for r in ctrl])
    cov = np.array([int(r[3]) for r in ctrl])
    assert 1700 < cov.min() and cov.max() < 2300, (cov.min(), cov.max())   # half B
    assert abs(frac.mean() - 0.02) < 0.004, frac.mean()
    assert [r[0] for r in ctrl] == ["oligoA"] * 112 + ["oligoB"] * 82
    assert ctrl[0][1] == str(SHIFT + 1) and ctrl[0][2] == GENOME["oligoA"][SHIFT]

    # Treated: modified positions recover their rates (a 10 SD shift is always
    # called; unshifted reads contribute the 2% background).
    trt = {(r[0], int(r[1])): r for r in out[("trt", "for")]}
    assert abs(float(trt[("oligoA", 41)][4]) - (0.60 + 0.40 * 0.02)) < 0.03
    assert abs(float(trt[("oligoA", 71)][4]) - (0.25 + 0.75 * 0.02)) < 0.03
    assert float(trt[("oligoA", 41)][6]) > 4 and float(trt[("oligoA", 71)][6]) < -1.5
    rev = np.array([float(r[4]) for r in out[("trt", "rev")]])
    assert rev.max() < 0.04, rev.max()                 # nothing on the other strand

    # Per-read table agrees with the summary.
    reads = read_table(os.path.join(tmp, "trt_for_reads.tsv.gz"))[1:]
    at41 = [r for r in reads if r[1] == "oligoA" and r[2] == "41"]
    assert len(at41) == int(trt[("oligoA", 41)][3])
    assert sum(float(r[6]) < 0.02 for r in at41) == int(trt[("oligoA", 41)][5])

    # A model from the other backend must be refused.
    try:
        run(perbase_mod_calls,
            ["--backend", "remora", "-i", bam, "-g", fasta, "-m", model, "-s", "for",
             "-o", os.path.join(tmp, "x.txt.gz")], lambda s, bam: treated[s])
    except SystemExit as exc:
        assert "backend" in str(exc)
    else:
        raise AssertionError("a model from another backend must be refused")

    # False-positive-rate threshold. With a Gaussian null it must land on the
    # p-value threshold: fpr 0.02 <-> |z| > 2.33, and the calls barely move.
    hists = []
    for strand in ("for", "rev"):
        hists.append(os.path.join(tmp, f"null_{strand}.hist.tsv"))
        run(perbase_mod_calls,
            common + ["-m", model, "-s", strand, "--half", "B", "--hist-out", hists[-1]],
            lambda s, bam: control[s])
    cutoff, n_null = mc.cutoff_from_hists(
        hists, 0.02, {"statistic": "absz", "fisher_lag": 0, "model": model})
    assert abs(cutoff - 2.326) < 0.03, cutoff
    assert 700000 < n_null < 850000, n_null
    fpr_path = os.path.join(tmp, "trt_for_fpr.txt.gz")
    run(perbase_mod_calls,
        common + ["-m", model, "-s", "for", "--fpr", "0.02", "--null-hist"] + hists
        + ["-o", fpr_path], lambda s, bam: treated[s])
    fpr = {(r[0], int(r[1])): r for r in read_table(fpr_path)}
    assert abs(float(fpr[("oligoA", 41)][4]) - float(trt[("oligoA", 41)][4])) < 0.005
    # A heavier-tailed null: the same fpr now needs a higher cutoff than 2.33,
    # which is the whole point of estimating it instead of assuming a Gaussian.
    heavy = np.zeros(mc.HIST_BINS + 1, dtype=np.int64)
    mc.hist_add(heavy, np.abs(np.random.default_rng(5).standard_t(3, 200000)))
    heavy_path = os.path.join(tmp, "heavy.hist.tsv")
    mc.write_hist(heavy_path, heavy, {"statistic": "absz", "fisher_lag": 0})
    assert mc.cutoff_from_hists([heavy_path], 0.02, {})[0] > 4.0
    # ... and a histogram from a different model is refused.
    try:
        mc.cutoff_from_hists(hists, 0.02, {"model": "other.tsv"})
    except SystemExit:
        pass
    else:
        raise AssertionError("a null histogram from another model must be refused")

    # Pooled model: two control BAMs whose levels differ by a constant offset (a
    # batch effect). The pooled mean sits between them and the pooled SD widens.
    def shifted(reads, offset):
        return [mc.ReadLevels(r.ref, r.start, r.read_id + "s", r.pos,
                              r.level + offset, r.kmer, r.clean) for r in reads]
    batch = {"a.bam": control["for"], "b.bam": shifted(control["for"], 0.2)}
    for name in batch:
        open(os.path.join(tmp, name), "w").close()
    pooled = os.path.join(tmp, "pooled.tsv")
    run(kmer_signal_model,
        ["--backend", "uncalled4", "-g", fasta, "--half", "A", "-o", pooled, "-i",
         os.path.join(tmp, "a.bam"), os.path.join(tmp, "b.bam")],
        lambda strand, bam: batch[os.path.basename(bam)])
    _, p_kmer, p_n, p_mean, p_sd = mc.load_model(pooled)
    assert p_kmer.size == 194
    assert np.allclose(p_mean, LEVEL[p_kmer] + 0.1, atol=0.03)
    assert np.allclose(p_sd, math.sqrt(NOISE ** 2 + 0.1 ** 2), rtol=0.15)

    # A sample with no reads on a reference (a B-DNA library on the PolyT reverse
    # complement) gives an empty model, and everything scored against it comes out
    # empty - without failing the workflow.
    empty_model = os.path.join(tmp, "empty_model.tsv")
    run(kmer_signal_model, common + ["--half", "A", "-o", empty_model],
        lambda strand, bam: [])
    assert mc.load_model(empty_model)[1].size == 0
    empty_calls = os.path.join(tmp, "empty_calls.txt.gz")
    empty_hist = os.path.join(tmp, "empty.hist.tsv")
    run(perbase_mod_calls,
        common + ["-m", empty_model, "-s", "for", "-o", empty_calls,
                  "--hist-out", empty_hist], lambda s, bam: treated[s][:50])
    assert read_table(empty_calls) == []
    cutoff, n_null = mc.cutoff_from_hists(
        [empty_hist], 0.02, {"statistic": "absz", "fisher_lag": 0})
    assert cutoff == float("inf") and n_null == 0
    # ... but reads that yield no clean k-mer at all are still an error.
    dirty = [mc.ReadLevels(r.ref, r.start, r.read_id, r.pos, r.level, r.kmer,
                           np.zeros(r.pos.size, dtype=bool))
             for r in control["for"][:20]]
    try:
        run(kmer_signal_model, common + ["--half", "all", "-o", empty_model],
            lambda strand, bam: dirty)
    except SystemExit as exc:
        assert "none had an error-free" in str(exc)
    else:
        raise AssertionError("reads with no clean k-mer must be an error")

    # Reactivity: merge across two contigs, treatment minus control.
    react = os.path.join(tmp, "react.txt.gz")
    old_argv = sys.argv
    sys.argv = ["mod_reactivity", "-p", os.path.join(tmp, "trt_for.txt.gz"),
                "-m", os.path.join(tmp, "ctrl_for.txt.gz"), "-g", fasta + ".fai",
                "-o", react]
    try:
        mod_reactivity.main()
    finally:
        sys.argv = old_argv
    rows = {(r[0], int(r[1])): r for r in read_table(react)}
    assert len(rows) == 194
    assert abs(float(rows[("oligoA", 41)][3]) - 0.588) < 0.04
    others = [float(r[3]) for key, r in rows.items()
              if key not in (("oligoA", 41), ("oligoA", 71))]
    assert max(abs(x) for x in others) < 0.03, max(abs(x) for x in others)


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and isinstance(fn, types.FunctionType):
            fn()
            print(f"ok  {name}")
    print("All tests passed.")
