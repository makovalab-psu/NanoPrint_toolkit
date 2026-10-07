#!/bin/bash

# run_case_tests.sh - Regression harness for the BAM-direct perbase_signal_deviation.
#
# For every edge-case BAM in cases/, builds the oracle the old way
# (Uncalled4_convert_tsv.sh -> the pre-rewrite perbase_signal_deviation.py), runs the
# new BAM-direct script, and compares. Both strands for each case.
#
# Two-tier verdict. First a per-column tolerance check with bounds derived from the
# oracle's %.6g rounding: 5e-6 for the mean and quantiles, more for mean_sq because
# d(x^2) = 2x*dx amplifies it by 2|x| and these BAMs are 1-6 reads deep, so nothing
# averages it away. If that is exceeded, the exact test runs before anything is called
# a failure: re-decode with %.6g applied and check byte-identity with the oracle. A
# case that fails the tolerance but proves byte-identical is a tolerance problem, and
# the harness says so rather than quietly passing or quietly failing.
#
# The oracle is built per strand because the old path required it: convert has no
# strand filter, so Uncalled4_convert_tsv.sh pre-filters the BAM with samtools. The
# new script does it in memory with -s.
#
# Needs: uncalled4, samtools, python3 + numpy + pysam.
#
# The test DATA is not in this repository -- it is ~104 MB of real mtDNA alignments
# and nine derived edge-case BAMs. Point -D at the directory holding it; see
# README.md in this folder for where it lives and how it was built.
#
# Usage, from anywhere:
#     workflow/tests/run_case_tests.sh -D <test_data_dir>

set -euo pipefail

# Toolkit root, derived from this script's own location (workflow/tests/ -> ../..),
# so the script works from any working directory.
TOOLKIT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
# Everything below is resolved against -D unless given as an absolute path.
DATA_DIR=""
CASE_DIR="cases"
GENOME="hg002v1.1_chrM.fa"
# The pre-rewrite script, which is the oracle generator. Pinned to the commit rather
# than HEAD: once the rewrite is committed, HEAD is the new script and the comparison
# would be against itself.
OLD_REF="7e72ff0"

usage() {
    cat << EOF
Usage: $(basename "$0") -D <test_data_dir> [-d <case_dir>] [-g <genome.fa>]
                           [-t <toolkit_dir>] [-r <git_ref>]

    -D  Directory holding the test data (required unless -d and -g are absolute).
        Expects <dir>/cases/ and <dir>/hg002v1.1_chrM.fa
    -d  Case BAM directory, relative to -D unless absolute (default: cases)
    -g  Reference FASTA, relative to -D unless absolute (default: hg002v1.1_chrM.fa)
    -t  NanoPrint_toolkit checkout (default: this script's own ../..)
    -r  Git ref holding the pre-rewrite script (default: ${OLD_REF})
    -h  This message
EOF
    exit 1
}

while getopts "D:d:g:t:r:h" opt; do
    case $opt in
        D) DATA_DIR="$OPTARG" ;;
        d) CASE_DIR="$OPTARG" ;;
        g) GENOME="$OPTARG" ;;
        t) TOOLKIT="$OPTARG" ;;
        r) OLD_REF="$OPTARG" ;;
        *) usage ;;
    esac
done

# Resolve the data paths against -D. Absolute paths win, so -d/-g can point
# anywhere without -D.
if [[ -n "$DATA_DIR" ]]; then
    [[ "$CASE_DIR" == /* ]] || CASE_DIR="${DATA_DIR%/}/${CASE_DIR}"
    [[ "$GENOME"   == /* ]] || GENOME="${DATA_DIR%/}/${GENOME}"
elif [[ "$CASE_DIR" != /* || "$GENOME" != /* ]]; then
    echo "Error: -D is required (or give -d and -g as absolute paths)." >&2
    echo "       The test data is not in this repository; see README.md here." >&2
    usage
fi

NEW_SCRIPT="${TOOLKIT}/workflow/scripts/perbase_signal_deviation.py"
CONVERT_SCRIPT="${TOOLKIT}/workflow/scripts/Uncalled4_convert_tsv.sh"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

for f in "$GENOME" "$NEW_SCRIPT" "$CONVERT_SCRIPT"; do
    [[ -e "$f" ]] || { echo "Error: not found: $f" >&2; exit 1; }
done
[[ -d "$CASE_DIR" ]] || { echo "Error: no case directory: $CASE_DIR" >&2; exit 1; }

command -v uncalled4 >/dev/null || { echo "Error: uncalled4 not on PATH" >&2; exit 1; }
command -v samtools  >/dev/null || { echo "Error: samtools not on PATH" >&2; exit 1; }

WORK="${CASE_DIR}/_oracle"
mkdir -p "$WORK"

# Recover the pre-rewrite script from git. If it is not in this checkout's history,
# say so plainly rather than silently comparing the new script against itself.
OLD_SCRIPT="${WORK}/old_perbase_signal_deviation.py"
if ! git -C "$TOOLKIT" show "${OLD_REF}:workflow/scripts/perbase_signal_deviation.py" \
        > "$OLD_SCRIPT" 2>/dev/null; then
    echo "Error: could not read workflow/scripts/perbase_signal_deviation.py at ${OLD_REF}" >&2
    echo "       Pass -r with a ref that predates the BAM-direct rewrite." >&2
    exit 1
fi
if grep -q 'direct from Uncalled4 BAM' "$OLD_SCRIPT"; then
    echo "Error: ${OLD_REF} already contains the rewritten script — it cannot be the" >&2
    echo "       oracle. Pass -r with an earlier ref." >&2
    exit 1
fi
echo "Oracle generator: ${OLD_REF} ($(wc -l < "$OLD_SCRIPT") lines)"
echo

PASS=0
FAIL=0
SKIP=0

for bam in "$CASE_DIR"/case*.bam; do
    name=$(basename "$bam" .bam)
    for strand in for rev; do
        label="${name} [${strand}]"
        tsv="${WORK}/${name}_${strand}.tsv"
        oracle="${WORK}/${name}_${strand}_oracle.txt.gz"
        new="${WORK}/${name}_${strand}_new.txt.gz"

        # An empty strand is legitimate for several cases (case6b is both-forward by
        # design) and is not a failure.
        n=$(samtools view -c $([[ "$strand" == "for" ]] && echo "-F" || echo "-f") 0x10 "$bam")
        if [[ "$n" -eq 0 ]]; then
            printf '  %-40s SKIP (no reads on this strand)\n' "$label"
            SKIP=$((SKIP + 1))
            continue
        fi

        "$CONVERT_SCRIPT" -i "$bam" -o "$tsv" -s "$strand" -g "$GENOME" -t 2 \
            > "${WORK}/${name}_${strand}.convert.log" 2>&1
        python3 "$OLD_SCRIPT" -i "$tsv" -g "$GENOME" -o "$oracle" \
            > "${WORK}/${name}_${strand}.oracle.log" 2>&1
        python3 "$NEW_SCRIPT" -i "$bam" -g "$GENOME" -s "$strand" -o "$new" \
            > "${WORK}/${name}_${strand}.new.log" 2>&1

        # Verdict: per-column agreement within the derived bounds. 5e-6 on the mean
        # and quantiles is the ceiling for a 6-significant-digit relative rounding;
        # mean_sq gets more room because d(x^2) = 2x*dx amplifies it by 2|x|, and these
        # BAMs are 1-6 reads deep so nothing averages it away.
        if python3 "${HERE}/compare_to_oracle.py" -a "$new" -b "$oracle" \
                > "${WORK}/${name}_${strand}.compare.log" 2>&1; then
            rows=$(gunzip -c "$new" | wc -l | tr -d ' ')
            sq=$(awk '/^mean_sq/ {print $4}' "${WORK}/${name}_${strand}.compare.log")
            printf '  %-40s PASS  rows=%-6s max mean_sq diff %s\n' \
                "$label" "$rows" "${sq:-n/a}"
            PASS=$((PASS + 1))
        else
            # Over tolerance. Escalate to the exact test before calling it a failure:
            # re-decode with the oracle's own %.6g quantisation and check byte-identity.
            # If that holds, the difference really is TSV rounding and the tolerance is
            # what is wrong -- which is worth saying out loud rather than hiding behind
            # a widened constant.
            if python3 "${HERE}/compare_to_oracle.py" --prove-precision \
                    -i "$bam" -g "$GENOME" -s "$strand" -b "$oracle" \
                    > "${WORK}/${name}_${strand}.prove.log" 2>&1; then
                printf '  %-40s PASS  over tolerance, but byte-identity PROVEN after\n' "$label"
                printf '  %-40s       %%.6g requantisation — widen the tolerance, not the code\n' ""
                PASS=$((PASS + 1))
            else
                printf '  %-40s FAIL  over tolerance AND not explained by TSV rounding\n' "$label"
                sed -n '/^column/,/^$/p' "${WORK}/${name}_${strand}.compare.log" | sed 's/^/      /'
                sed -n '1,6p' "${WORK}/${name}_${strand}.prove.log" | sed 's/^/      /'
                FAIL=$((FAIL + 1))
            fi
        fi

        rm -f "$tsv"        # the TSV is the big file; the outputs are small
    done
done

echo
echo "PASS: $PASS   FAIL: $FAIL   SKIP: $SKIP"
echo "Logs and outputs kept in ${WORK}/"
[[ "$FAIL" -eq 0 ]]
