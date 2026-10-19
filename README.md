# C_{3a} > 1.19281555: a finite k = 2 controlled-carry instance (DRAFT, local only, not yet published)

**Author / responsible human:** Sylla ([@xSyllaa](https://github.com/xSyllaa)).
**Status:** draft, not yet published or submitted. On 2026-10-19 three adversarial replays run inside the same AI pipeline (fresh contexts, each with its own code and
its own histogram rebuild) confirmed the value, and a hostile review found no mathematical gap. This is **not** an
independent third-party replay in the sense of CONTRIBUTING; there has been no external review yet.
Hosting (public repo and/or Zenodo DOI + pinned commit) still to be decided by Sylla.

## Result
For the Gyarmati–Hennecart–Ruzsa sum–difference constant ([page 3a](https://teorth.github.io/optimizationproblems/constants/3a.html)),

  **C_{3a} > 1.19281555**

from one explicitly specified finite set V = 2U (2×10^18 base-27022 digits), via the GHR lemma.
This is a **finite instantiation of Kleinwaks' controlled-carry block bound** [K2026b, Prop. 3] at block length k = 2,
using K2026b's record mask ⟨1971,2016,2100,2628,2688,2800⟩∩[0,26972] and base 27022, with three explicit choices:
x = 0x1.ffde1bafccc06p-1 (optimised for k = 2), ℓ = 10^18 blocks, and no blocks at the self-inverse residue q²/2 (N_h = 0).
It is not a new construction; the idea and the asymptotic proposition are Kleinwaks'. Open PR #236 (mizopl) is an earlier
finite instantiation of the controlled-carry idea (k = 1, odd base, 1.19230254).

Limit of the same family at this x (Comments column, asterisked, not claimed): 1.19281555778…*.
1.19281555 is the largest 8-decimal value this certificate chain can give at k = 2 with this mask and base (proof/THEOREM.md §5).

## Verify (≈1 minute)
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
| `SHA256SUMS` | hashes of every other file in this repository |
| `LICENSE` | MIT |

Exact data: n₀ = 3362808316958, C = 26096544438312721237750, |M| = 5869, |M+M| = 31274, θ_lo = 1.192815557092986…

## AI-use disclosure and responsibility
This work was carried out with an AI agent pipeline (Grok Bot) directed by Sylla: it searched parameters, wrote and ran
the checkers, drafted this package and the theorem note, and ran adversarial "skeptic" reviews of the claim.
Sylla (xSyllaa) directed the work, reviewed the argument and the outputs, and takes responsibility for the result.
All references were checked against the cited sources (K2026b at commit 352ab86, the 3a page and CONTRIBUTING.md).

## References
- [K2026b] L. Kleinwaks, *Improved lower bound for the Gyarmati–Hennecart–Ruzsa sum–difference constant using masked digits and controlled carries*,
  github.com/kleinwaks/masked-digit-sum-difference-bound (commit 352ab86), Prop. 3 "Controlled-carry block bound".
- [GHR2007] K. Gyarmati, F. Hennecart, I. Z. Ruzsa, *Sums and differences of finite sets*, Funct. Approx. 37 (2007) 175–186.
- H. Robbins, *A remark on Stirling's formula*, Amer. Math. Monthly 62 (1955) 26–29.
- PR #236 (mizopl), teorth/optimizationproblems, finite controlled-carry instance at k = 1.

## Credit
Sylla / xSyllaa. Checker lineage includes code written during adversarial reviews in this pipeline.
