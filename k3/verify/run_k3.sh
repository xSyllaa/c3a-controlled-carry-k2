#!/usr/bin/env bash
# k=3 replay. Quick (default, ~1 min): checks the shipped k=3 histogram and runs the certificates.
# FULL=1: also rebuilds the histogram twice from source (~20 min on 8 threads, ~6 GB RAM, gcc with OpenMP)
# and checks that builder A's k=2 by-product equals the published k=2 histogram.
# Needs gcc, python3, libmpfr.so.6; optional mpmath + numpy for the secondary checker. Credit: Sylla / xSyllaa.
set -euo pipefail
cd "$(dirname "$0")/.."
ELL=1000000000000000000000000; X=0x1.ffde3d0d5025dp-1; CLAIM=1.19359680
EXP3=3186f6e220fd4a36369bf23ce841a9d54e6179cfaff00b5f3189c058f0d85c63
EXP2=4d16010f1b41d875c2912f397e76e8d4010e49b635e139ad1f9bcbd741956c63
TMP=$(mktemp -d); H="$PWD/data/k3_cost_hist.tsv"
[ "$(sha256sum $H | cut -d' ' -f1)" = "$EXP3" ] && echo "sha256 ok  $H" || { echo "FAIL sha $H"; exit 1; }
if [ "${FULL:-0}" = 1 ]; then
  echo "== rebuild A (carry automaton, low-heavy)"; gcc -O3 -fopenmp -o "$TMP/A" src/k3_auto_A.c -lm; "$TMP/A" "$TMP/k3A.tsv" "$TMP/k2A.tsv" 2>/dev/null
  echo "== rebuild B (c2 by DxD scatter, top-heavy)"; gcc -O3 -fopenmp -o "$TMP/B" src/k3_pairs_B.c -lm; "$TMP/B" "$TMP/k3B.tsv" 2>/dev/null
  for f in "$TMP/k3A.tsv" "$TMP/k3B.tsv"; do [ "$(sha256sum "$f" | cut -d' ' -f1)" = "$EXP3" ] || { echo "FAIL sha $f"; exit 1; }; echo "sha256 ok  $f"; done
  [ "$(sha256sum "$TMP/k2A.tsv" | cut -d' ' -f1)" = "$EXP2" ] || { echo "FAIL k2 by-product"; exit 1; }; echo "sha256 ok  k=2 by-product of A = ../data/k2_cost_hist.tsv"
  H="$TMP/k3A.tsv"
fi
echo "== PRIMARY k=3 certificate (MPFR directed rounding + Robbins)"
(cd verify && python3 cert_k3_primary.py --hist "$H" --json-out "$TMP/p3.json" | tail -1)
echo "== regression: same script at k=2 reproduces the PR #238 certificate"
(cd verify && python3 cert_k3_primary.py --k 2 --hist ../../data/k2_cost_hist.tsv --ell 1000000000000000000 --xhex 0x1.ffde1bafccc06p-1 --claim 1.19281555 --json-out "$TMP/p2.json" | tail -1)
python3 - "$TMP/p2.json" ../certificates/cert_primary_ell_1e18.json <<'PY'
import json,sys
a,b=(json.load(open(p)) for p in sys.argv[1:3])
assert a["n0"]==str(b["n0"]) and a["C"]==str(b["C"]) and a["theta_lo"][:40]==str(b["theta_lo"])[:40], "k=2 regression mismatch"
print("k=2 regression ok (n0, C, theta_lo identical to certificates/cert_primary_ell_1e18.json)")
PY
echo "== negative control (claim 1.19359681 must fail)"
if (cd verify && python3 cert_k3_primary.py --hist "$H" --claim 1.19359681 >/dev/null); then echo "FAIL negative control passed"; exit 1; else echo "negative control ok"; fi
if python3 -c 'import mpmath, numpy' 2>/dev/null; then
  echo "== secondary (mpmath.iv intervals, exact low-cost prefix, sampled residues)"
  (cd verify && python3 cert_k3_secondary_mpmath_iv.py "$H" $ELL $X $CLAIM | grep -q '"PASS": true' && echo "secondary PASS")
  SEC=run
else SEC="skipped (pip install mpmath numpy)"; fi
echo "ALL k=3 CHECKS PASS: C_3a > $CLAIM   (secondary: ${SEC})"
