# C_{3a} > 1.19359680 (finite k = 3) and C_{3a} > 1.19281555 (finite k = 2): controlled-carry instances

**Author / responsible human:** Sylla ([@xSyllaa](https://github.com/xSyllaa)).
**Status:** public verification package for teorth/optimizationproblems
[PR #238](https://github.com/teorth/optimizationproblems/pull/238). The k = 2 instance (1.19281555) was submitted on
2026-10-19 at commit ad16a46; the k = 3 instance (1.19359680, found 2026-10-20) is in `k3/`. The k = 2 files are unchanged
and still reproduce on their own. Not peer-reviewed and not externally reviewed. For each instance, three adversarial
replays run inside the same AI pipeline (separate processes with their own code and histogram rebuilds) confirmed the
value, and a hostile review found no mathematical gap. This is **not** an independent third-party replay in the sense of
CONTRIBUTING. Hosted at https://github.com/xSyllaa/c3a-controlled-carry-k2 (the repository name predates the k = 3 instance; no Zenodo DOI).

## Result
For the Gyarmati–Hennecart–Ruzsa sum–difference constant ([page 3a](https://teorth.github.io/optimizationproblems/constants/3a.html)),

  **C_{3a} > 1.19359680**   (k = 3: one explicit finite set V = 2U with 3×10^24 base-27022 digits; `k3/proof/THEOREM_k.md`)
  C_{3a} > 1.19281555   (k = 2: 2×10^18 digits; `proof/THEOREM.md`)

Both are **finite instantiations of Kleinwaks' controlled-carry block bound** [K2026b, Prop. 3], at block lengths k = 3
and k = 2, with K2026b's record mask ⟨1971,2016,2100,2628,2688,2800⟩∩[0,26972] and base 27022. Our choices are x
(optimised per k), ℓ (10^24 resp. 10^18 blocks), and no blocks at the self-inverse residue q^k/2 (N_h = 0). It is not a
new construction; the idea and the asymptotic proposition are Kleinwaks'. Open PR #236 (mizopl) is an earlier finite
instantiation of the controlled-carry idea (k = 1, odd base, 1.19230254).

The new ingredient at k = 3 is the exact min-cost histogram mod q³ (q³ ≈ 1.97×10^13 residues), computed by a 2-state
carry automaton (Lemma 0 of `k3/proof/THEOREM_k.md` proves it gives the true minimum) factorised into joint cost pairs.

Limits of the same family (Comments column, asterisked, not claimed): k = 3: 1.19359680913…*; k = 2: 1.19281555778…*.
1.19359680 is the largest 8-decimal value this certificate chain gives at k = 3 with this mask and base.

## Verify k = 3 (≈1 minute; full rebuild ≈20 minutes, ~6 GB RAM)
```bash
sudo apt install build-essential libmpfr6
python3 -m pip install mpmath numpy              # optional, for the secondary checker
bash k3/verify/run_k3.sh                         # last line: ALL k=3 CHECKS PASS: C_3a > 1.19359680 (secondary: run)
FULL=1 bash k3/verify/run_k3.sh                  # also rebuilds the k=3 histogram twice from source
```
Our outputs are in `k3/verify/run_k3_output.txt` and `k3/verify/run_k3_full_output.txt`.

## Verify k = 2 (≈1 minute)
```bash
sudo apt install build-essential libmpfr6        # libgmp comes with it
python3 -m pip install mpmath                    # optional, for the two secondary checkers
bash verify/run_all.sh                           # last line: ALL CHECKS PASS: C_3a > 1.19281555 (secondary checkers: run)
```
Our output is in `verify/run_all_output.txt`. Requires Linux (the MPFR binding loads `libmpfr.so.6`). Primary checker alone (stdlib python + libmpfr, ~5 s):
`cd verify && python3 cert_primary.py --ell 1000000000000000000 --xhex 0x1.ffde1bafccc06p-1 --claim 1.19281555`.

## Contents
| path | what |
|---|---|
| `proof/THEOREM.md` | the finite theorem and its lemma chain, with every explicit number |
| `src/finite_k2.c` | k=2 min-cost residue histogram (full q² array, ~1.5 GB) |
| `src/hist_indep.c` | independent histogram builder (different lift enumeration, light memory) |
| `data/k2_cost_hist.tsv` | 26540-row cost histogram, sha256 `4d16010f1b41d875c2912f397e76e8d4010e49b635e139ad1f9bcbd741956c63` |
| `data/sumset.txt` | M+M (31274 values) |
| `data/k2_meta.json` | generator metadata: x-independent facts (B, q, Q2, c(0), c(h), 26540 costs, max cost 53762) and the certificate x; its Z2/P+ fields are binary64 info only (certified enclosures are in `certificates/`) |
| `verify/cert_primary.py` + `mpfr_ctypes.py` | **primary**: MPFR directed rounding, Robbins ln n! bounds, type counts from enclosures |
| `verify/cert_secondary_mpmath.py` | secondary: mpmath 120-digit loggamma |
| `verify/cert_indep_robbins.py` | independent: mpmath Robbins on the independently rebuilt histogram |
| `verify/toy_bruteforce.py` | exact enumeration of the same construction for q = 4, 6 |
| `verify/xopt_solve.py` | x-stationarity solve (optimisation only; any 0<x<1 is valid) |
| `certificates/*.json` | outputs at ℓ = 10^18 (and 10^20 cross-check) |
| `k3/proof/THEOREM_k.md` | block-length-k theorem note (k = 2 and k = 3 instances), incl. Lemma 0 (carry automaton = true min) |
| `k3/src/k3_auto_A.c` | k=3 histogram builder A: low-heavy carry automaton (also emits the k=2 histogram, = `data/k2_cost_hist.tsv`) |
| `k3/src/k3_pairs_B.c` | k=3 histogram builder B: top-heavy, c₂ by full D×D scatter |
| `k3/data/k3_cost_hist.tsv` | 54099-row k=3 cost histogram, sha256 `3186f6e220fd4a36369bf23ce841a9d54e6179cfaff00b5f3189c058f0d85c63` |
| `k3/verify/cert_k3_primary.py` + `mpfr_ctypes.py` | **primary k=3**: MPFR directed rounding + Robbins; `--k 2` reproduces the k=2 certificate |
| `k3/verify/cert_k3_secondary_mpmath_iv.py` | secondary: mpmath interval arithmetic, exact low-cost prefix, sampled residues |
| `k3/verify/run_k3.sh` | k=3 replay (quick / `FULL=1`) |
| `k3/certificates/cert_k3_primary_ell_1e24.json` | k=3 certificate output |
| `SHA256SUMS` | hashes of every other file in this repository |
| `LICENSE` | MIT |

Exact data (k = 3): n₀ = 5794476231064238, C = 39290318538006294714264853344, 54099 distinct costs, θ_lo = 1.193596809129730…
Exact data (k = 2): n₀ = 3362808316958, C = 26096544438312721237750, |M| = 5869, |M+M| = 31274, θ_lo = 1.192815557092986…

## AI-use disclosure and responsibility
This work was carried out with an AI agent pipeline (Grok Bot) directed by Sylla: it searched parameters, wrote and ran
the checkers, drafted this package and the theorem note, and ran adversarial "skeptic" reviews of the claim.
Sylla (xSyllaa) directed the work and takes responsibility for the result.
The pipeline checked all references against the cited sources (K2026b at commit 352ab86, the 3a page and CONTRIBUTING.md).

## References
- [K2026b] L. Kleinwaks, *Improved lower bound for the Gyarmati–Hennecart–Ruzsa sum–difference constant using masked digits and controlled carries*,
  github.com/kleinwaks/masked-digit-sum-difference-bound (commit 352ab86), Prop. 3 "Controlled-carry block bound".
- [GHR2007] K. Gyarmati, F. Hennecart, I. Z. Ruzsa, *Sums and differences of finite sets*, Funct. Approx. 37 (2007) 175–186.
- H. Robbins, *A remark on Stirling's formula*, Amer. Math. Monthly 62 (1955) 26–29.
- PR #236 (mizopl), teorth/optimizationproblems, finite controlled-carry instance at k = 1.


## Credit
Sylla / xSyllaa. Checker lineage includes code written during adversarial reviews in this pipeline.
