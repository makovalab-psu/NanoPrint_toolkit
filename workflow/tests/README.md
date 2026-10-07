# Regression tests for `perbase_signal_deviation`

These test the DTW signal branch — `perbase_signal_deviation.py` reading an Uncalled4
BAM directly. The oracle is the path it replaced: `uncalled4 convert` → the pre-rewrite
script. Keeping `uncalled4_convert_tsv` in the workflow is partly what makes that
possible.

| Script | What it does |
|---|---|
| `make_dtw_test_cases.py` | Derives nine edge-case BAMs from a real Uncalled4 BAM |
| `compare_to_oracle.py` | Per-column agreement against an oracle output |
| `run_case_tests.sh` | Drives every case through convert + the old script + the comparison |

## The test data is not in this repository

It is ~104 MB: a real 3,171-read chrM Uncalled4 BAM, its `perbase_signal` output as the
oracle, the chrM reference, and the nine derived case BAMs. It lives with the project
that produced it:

```
js4022_Permanganate_probing_protocol_in_bacteria/dtw_signal_calculation_test_data/
```

Point `-D` at that directory. Its own `README.md` records the provenance in detail.

**Where the BAM came from, and the trap in getting it.** It is the full-depth chrM
extraction from js4028:

```
js4028_Subsample_mtDNA_reads_and_recalculate_error_rates_2/NanoPrint_toolkit/data/mtDNA/
```

Copy from `data/mtDNA/`, **not** `data/uncalled4/` — in js4028 the latter holds only
staging symlinks into `../../mtDNA/` (`01_subsample_mtDNA.sbatch:394-399`, the pre-`^u`
workaround). A copy that does not dereference them yields a truncated file: the first
attempt produced a 42-byte BAM whose BGZF header declared a 968-byte block. Use
`rsync -L` and `samtools quickcheck` the result.

## Running

```bash
conda activate <env with uncalled4, samtools, numpy, pysam>
workflow/tests/run_case_tests.sh -D /path/to/dtw_signal_calculation_test_data
```

Expect a PASS line per case per strand and `FAIL: 0`. SKIP is legitimate where a case
has no reads on that strand — `case6b` is deliberately both-forward, because an all-NA
position has to be all-NA *within* a strand to test anything.

The harness recovers the oracle generator from git (`-r`, default `7e72ff0`, the last
commit before the BAM-direct rewrite) and **refuses to run if that ref already contains
the rewrite** — otherwise it would compare the new script against itself and pass for
the wrong reason.

## What the cases cover

The mtDNA oracle exercises the main path well: both strands, 16,560 positions each,
~250–666× depth. It does not reach the awkward shapes, which is what `cases/` is for.
All nine occur naturally in the real BAM — none are synthesised from scratch, so
uncalled4's k-mer registration comes along for free instead of being re-derived (and
therefore possibly re-derived *wrongly* in both the fixture and the code under test).

| Case | Covers |
|---|---|
| `case1_both_strands` | fwd+rev over the same positions — the `mpos` negation and swapped k-mer shift |
| `case2_multi_interval` | `ur` with several intervals (reference deletion) |
| `case3_na_positions` | NA sentinel in `uc`, one read per strand |
| `case4_skips` | `ul == 0`, consecutive positions sharing a sample interval |
| `case5_window_boundary` | reads crossing a window boundary — no double-count, no drop |
| `case6a_coverage1` | coverage 1 within a strand (quantile edge) |
| `case6b_all_na_position` | a stretch blanked in every read covering it |
| `case7_split_read` | the `pi` tag, where `sam_to_aln` branches |
| `case8_long_runlen` | `<= MIN_I16` run-length markers |

## Expected disagreement with the oracle

The two paths are **not** byte-identical, and that is not a defect in either.
`uncalled4 convert` writes its TSV with `float_format="%.6g"` (`src/uncalled4/io/tsv.py`),
so the oracle's statistics come from values already rounded to six significant digits.
The BAM-direct path never serialises and keeps the full float64.

What must match **exactly**: row count, chromosome, position, nucleotide, coverage. Those
come from the decode rather than the arithmetic, and `compare_to_oracle.py` hard-fails on
any mismatch there.

What may differ: the mean, the four quantiles and `mean_sq`, in the last printed digits.
Two separate bounds, because they are not the same size:

| Columns | Bound | Why |
|---|---|---|
| mean, q25, q75, q025, q975 | `--max-abs`, default **1e-5** | `%.6g` is a relative error of ≤5e-6, so ~5e-6 absolute near 1. Observed ceiling on real data: exactly 5e-6 |
| mean_sq | `--max-abs-sq`, default **5e-5** | `d(x²) = 2x·dx`, so squaring amplifies each observation's rounding by **2\|x\|** — measured 2× to 20× for \|x\| from 0.5 to 2.5 |

**Coverage is the other half of it.** Averaging cancels these errors, so deep data agrees
far more closely than shallow:

| Data | Depth | Worst `mean_sq` difference |
|---|---|---|
| Full mtDNA sample | 250–666× | **1e-6** |
| The nine case BAMs | 1–6 reads | **2.5e-5** |

The cases are a deliberate worst case — one read means `mean_sq` *is* x², with nothing to
average the rounding away. **Do not read a larger spread on the cases as a regression.**
5e-5 is still only ~0.2% of a typical `mean_sq` (median 0.024 in real data), so it cannot
hide an actual arithmetic defect.

**If the tolerance is exceeded, the harness escalates rather than failing outright.** It
re-runs the decode with `%.6g` applied to every observation — exactly what
`uncalled4 convert` does — and checks byte-identity with the oracle. A case that is over
tolerance but proves byte-identical is a *tolerance* problem, and the harness says so
instead of quietly passing or quietly failing. Only a case that is both over tolerance
and not reproducible by requantisation is a real FAIL.

**Where they differ, the BAM-direct path is the more accurate one. Do not tune it toward
the oracle.**

## History

First full run, 2026-10-07: **all nine cases decoded correctly on both strands** —
`chr / pos / nt / coverage` identical on every row, which is what certifies the k-mer
registration, the reverse-strand `mpos` flip, the strand filter and the NaN handling
across deletions, skips, split reads, NA positions, window boundaries, coverage-1 and
long run-length markers.

That run reported 14 FAILs, all of them the tolerance being wrong rather than the code:
a single `--max-abs 1e-5` was applied to all six numeric columns, derived from the linear
bound, with no allowance for the `2|x|` amplification in `mean_sq`. Every observed value
was consistent with pure `%.6g` quantisation. Fixed by splitting the tolerance per column
and adding the escalation above.

**After the fix: `PASS: 17  FAIL: 0  SKIP: 1`.** Observed `mean_sq` spread 9e-6 to
2.5e-5, against the 5e-5 bound — at least 2x headroom on every case, and **no case
needed the requantisation escalation**, so the bound is calibrated rather than merely
loose. Worst cases are `case2_multi_interval [for]` at 2.5e-5 and
`case1_both_strands [rev]` at 2.3e-5; treat anything materially above 2.5e-5 on this
data as worth investigating even though it would still pass.
