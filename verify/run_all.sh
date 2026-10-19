#!/usr/bin/env bash
# One-command replay. Needs: gcc, python3, libmpfr.so.6 + libgmp.so.10 (Debian/Ubuntu: apt install libmpfr6);
# optional: python3 -m pip install mpmath  (for the two secondary checkers).
# Credit: Sylla / xSyllaa.
set -euo pipefail
cd "$(dirname "$0")/.."
ELL=1000000000000000000; X=0x1.ffde1bafccc06p-1; CLAIM=1.19281555
TMP=$(mktemp -d)
echo "== 1. rebuild cost histogram from source (finite_k2.c, ~40 s, ~1.5 GB RAM)"
gcc -O2 -o "$TMP/finite_k2" src/finite_k2.c -lm
"$TMP/finite_k2" "$TMP/k2_cost_hist.tsv" "$TMP/meta.json" 2>/dev/null
python3 - "$TMP/meta.json" data/k2_meta.json <<'PY'
import json, sys
a, b = (json.load(open(p)) for p in sys.argv[1:3])
keys = ["B", "q", "k", "Q2", "covered", "c0", "ch", "distinct_costs", "max_cost", "x_hex_certificate"]
bad = [k for k in keys if a.get(k) != b.get(k)]
if bad: sys.exit(f"FAIL meta mismatch on {bad}")
print("meta ok    (x-independent fields + certificate x match data/k2_meta.json)")
PY
echo "== 2. independent histogram builder (hist_indep.c)"
gcc -O2 -o "$TMP/hist_indep" src/hist_indep.c
(cd "$TMP" && ./hist_indep 2>/dev/null)
EXP=4d16010f1b41d875c2912f397e76e8d4010e49b635e139ad1f9bcbd741956c63
for f in data/k2_cost_hist.tsv "$TMP/k2_cost_hist.tsv" "$TMP/my_hist.tsv"; do
  h=$(sha256sum "$f" | cut -d' ' -f1); [ "$h" = "$EXP" ] || { echo "FAIL sha $f $h"; exit 1; }; echo "sha256 ok  $f"; done
cmp -s data/sumset.txt "$TMP/sumset.txt" && echo "sumset ok  (data/sumset.txt = independently rebuilt M+M)" || { echo "FAIL sumset"; exit 1; }
echo "== 3. PRIMARY certificate (MPFR directed rounding + Robbins; stdlib python)"
(cd verify && python3 cert_primary.py --ell $ELL --xhex $X --claim $CLAIM --hist "$TMP/k2_cost_hist.tsv" --json-out "$TMP/primary.json" | tail -1)
echo "== 4. toy brute force of the construction (q=4,6)"
(cd verify && python3 toy_bruteforce.py | tail -1)
SEC="skipped (pip install mpmath to run them)"
if python3 -c 'import mpmath' 2>/dev/null; then
  echo "== 5. secondary checkers (mpmath)"
  (cd verify && python3 cert_secondary_mpmath.py --ell $ELL --xhex $X --claim $CLAIM >/dev/null && echo "secondary PASS")
  (cd verify && python3 cert_indep_robbins.py $ELL $X "$TMP/my_hist.tsv" "$TMP/sumset.txt" > "$TMP/indep.json")
  python3 - "$TMP/primary.json" "$TMP/indep.json" <<'PY'
import json, sys
p, i = (json.load(open(f)) for f in sys.argv[1:3])
assert i["gt_1_19281555"] is True, "independent checker: theta <= claim"
assert str(i["n0"]) == p["n0"] and str(i["C"]) == p["C"], "independent n0/C differ from primary"
assert i["floor8"] == p["theta_floor8"] == "1.19281555", "floor8 mismatch"
print(f'independent PASS  (n0={i["n0"]}, C={i["C"]}, floor8={i["floor8"]}, theta={i["theta"]})')
PY
  SEC="run"
else echo "== 5. skipped (pip install mpmath to run secondary checkers)"; fi
echo "ALL CHECKS PASS: C_3a > $CLAIM   (secondary checkers: $SEC)"
