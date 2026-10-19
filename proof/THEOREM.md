# Finite theorem note: C_{3a} > 1.19281555 (k = 2 controlled-carry, N_h = 0)

Author and responsible human: Sylla (xSyllaa). Checked on 2026-10-19 by three adversarial replays inside the
same AI pipeline (fresh contexts, own code, own histogram rebuilds) and one hostile review; not yet externally reviewed. This note states and proves the *finite* instance that
the certificate checks. It is a finite instantiation of Kleinwaks' controlled-carry block bound
[K2026b, Prop. 3 `prop:controlled-block`, repo kleinwaks/masked-digit-sum-difference-bound @352ab86]
at k = 2 with the K2026b record mask and base. The only departures from K2026b's (asymptotic)
proof are: a fixed explicit ℓ, an explicit x, **no h-blocks (N_h = 0)**, and L = C/2 exactly
(no ε). No new construction is claimed.

## 0. Data
- Generators G = {1971, 2016, 2100, 2628, 2688, 2800}; mask M = ⟨G⟩ ∩ [0, B], B = 26972.
  |M| = 5869, 0 ∈ M, B ∈ M, |M+M| = 31274.
- Base q = 27022 (> B), block length k = 2, Q = q² = 730188484. Self-inverse residues mod Q: 0 and h = Q/2.
- x = 0x1.ffde1bafccc06p-1 = 9004870224628742 / 2^53 ≈ 0.99974142571443 (exact dyadic, 0 < x < 1).
- ℓ = 10^18 blocks, so N = 2ℓ = 2×10^18 base-q digits.

## 1. Costs and lifts
For d ∈ D = M − M let κ(d) = min{a + b : a, b ∈ M, a − b = d}, with witnesses (α_d, β_d) attaining
it, chosen so that α_{−d} = β_d, β_{−d} = α_d, α_0 = β_0 = 0.

For a residue r ∈ Z/QZ let c(r) = min{κ(d₀) + κ(d₁) : d_i ∈ D, d₀ + d₁q ≡ r (mod Q)}.
Every residue is attainable (all Q residues are covered). The **cost histogram**
m(c) = #{r : c(r) = c} has 26540 distinct costs, min 0 (only r = 0), max 53762, c(h) = 24367.
It is fully determined by (G, B, q) and is recomputed by two independent programs
(`finite_k2.c`: full Q-array min over D×D; `hist_indep.c`: per-digit carry-split rule),
both giving `k2_cost_hist.tsv` with sha256 `4d16010f1b41d875c2912f397e76e8d4010e49b635e139ad1f9bcbd741956c63`.
Since c(−r) = c(r) (negate the digits), every cost class except those of 0 and h has even size.

**Lifts.** Fix one representative r of every orbit {r, −r} with r ∉ {0, h}, choose a min-cost
vector (d₀(r), d₁(r)) and put (d₀(−r), d₁(−r)) = −(d₀(r), d₁(r)); for 0 take (0, 0).
The integer value e(r) = d₀(r) + d₁(r) q satisfies e(r) ≡ r (mod Q), e(−r) = −e(r), and e is a
*function of the residue* (this is where N_h = 0 is used: h is never used, so no orientation choice arises).
Block realisations: u-block (α_{d₀}, α_{d₁}), v-block (β_{d₀}, β_{d₁}); their difference has value e(r);
swapping gives the realisation of −r.

## 2. Type and the set U
Z₂(x) = Σ_r x^{c(r)} = Σ_c m(c) x^c. Type counts (K2026b's rule, with h removed):
  N_r = N_{−r} = ⌊ℓ x^{c(r)} / Z₂(x)⌋ for r ∉ {0, h};  N_h = 0;  N₀ = ℓ − Σ_{r≠0,h} N_r.
(N_r depends only on c(r), so N_r = N_{−r}.) Exact values: **n₀ = N₀ = 3362808316958 ≥ 0**,
all 730188482 residues r ∉ {0, h} have N_r ≥ 1, **C = Σ_r N_r c(r) = 26096544438312721237750** (even),
**L = C/2 = 13048272219156360618875**.

U = U_ℓ = { Σ_{j<2ℓ} a_j q^j : a_j ∈ M, Σ_j a_j ≤ L } ⊂ Z_{≥0}, 0 ∈ U.

## 3. Lemmas
**Lemma A (encoding lands in U).** For a word w ∈ (Z/QZ)^ℓ of the type above put
u(w) = Σ_i (α_{d₀(w_i)} + α_{d₁(w_i)} q) Q^i, v(w) = Σ_i (β_{d₀(w_i)} + β_{d₁(w_i)} q) Q^i.
All base-q digits of u, v lie in M (each ≤ B < q, no carries). Since r and −r occur equally often and
the realisation of −r is the swap of that of r, digit-sum(u) = digit-sum(v) = Σ_{pairs} N_r c(r) = C/2 = L.
So u(w), v(w) ∈ U, and u(w) − v(w) = Σ_i e(w_i) Q^i.

**Lemma B (injectivity).** w ↦ Δ(w) = Σ_i e(w_i) Q^i is injective: Δ mod Q = w₀ determines w₀,
hence e(w₀) (e is a function of the residue); (Δ − e(w₀))/Q = Δ(w₁w₂…) and induct.
(|e(r)| ≤ B(q+1) = 728864356 < Q, as in K2026b; not even needed for this argument.)
Hence |U − U| ≥ #words = ℓ! / (N₀! Π_{r≠0,h} N_r!) =: 𝓜.

**Lemma C (sumset).** U + U ⊆ { Σ_{j<2ℓ} s_j q^j : s_j ∈ M+M, Σ s_j ≤ 2L = C }. For 0 < x < 1,
|U + U| ≤ #{s ∈ (M+M)^{2ℓ} : Σ s_j ≤ C} ≤ Σ_s x^{Σ s_j − C} = x^{−C} P₊(x)^{2ℓ}, P₊(x) = Σ_{s ∈ M+M} x^s.

**Lemma D (GHR with V = 2U).** [GHR2007] (as quoted on the 3a page): a finite set V ⊂ Z_{≥0} with 0 ∈ V and
|V − V| < 2maxV + 1 gives C_{3a} ≥ 1 + log(|V−V|/|V+V|) / log(2maxV + 1). Take V = 2U:
V − V = 2(U − U) has |V−V| = |U−U| ≤ 2maxU + 1 < 4maxU + 1 = 2maxV + 1 (strict since maxU > 0), and |V+V| = |U+U|.
maxU ≤ B(q^{2ℓ} − 1)/(q − 1), so log(2maxV+1) ≤ log(4B q^{2ℓ}/(q−1) + 1).

## 4. Theorem
With the data of §0–§2,
  C_{3a} ≥ 1 + (log 𝓜 − 2ℓ log P₊(x) + C log x) / log(4B q^{2ℓ}/(q−1) + 1)
(the numerator is positive, so replacing |U−U|, |U+U|, the denominator by the bounds of Lemmas B–D only
lowers the right side). Evaluated with MPFR directed rounding and Robbins' ln n! bounds (`cert_primary.py`):

| quantity | certified bound |
|---|---|
| log 𝓜 ≥ | 19353442719252255618.41914067403… |
| log(x^{−C} P₊^{2ℓ}) ≤ | 15418306021682217276.01171204408… |
| log(4Bq^{2ℓ}/(q−1)+1) ≤ | 20408813256039808818.33351550874… |
| Z₂(x) ≥ | 297948.172635433171794… (√Z₂ ≈ 545.84629030106) |
| P₊(x) ≤ | 76.30733448716651847590… |
| **θ ≥** | **1.192815557092986248985…** |

**Therefore C_{3a} > 1.19281555.** (θ_floor8 = 1.19281555, θ_floor12 = 1.192815557092.)

## 5. Limit (for the Comments column, asterisked, not claimed)
As ℓ → ∞ at this x the same family tends to 1 + (½ log Z₂(x) − log P₊(x)) / log q ∈
[1.1928155577886288831, …832]* (= K2026b Prop. 3 at k = 2). x is the stationary point of ½E_p[c] = E_{P₊}[s] (|g(x)| < 1e-9, loss vs exact optimum < 1e-27).
*Why no further (x, ℓ, count-rule) bump exists.* Use any Chernoff parameter x' ∈ (0,1) in Lemma C and any symmetric
N_h = 0 type p (frequencies N_r/ℓ). Then log 𝓜 ≤ ℓH(p) and C = ℓE_p[c], so the numerator of §4 is at most
ℓ[H(p) + E_p[c] log x' − 2 log P₊(x')] ≤ ℓ[log Z₂(x') − 2 log P₊(x')] (Gibbs: max_p [H(p) + E_p[c] log x'] = log Z₂(x')),
while the denominator is ≥ 2ℓ log q. So the right side of §4 is < θ∞(x') := 1 + (½ log Z₂(x') − log P₊(x'))/log q for every ℓ.
Numerically θ∞ is unimodal on (0,1) with its maximum at this x: a log-grid scan of t = −log x' over [1e-7, 100]
and a fine scan of ±50% around t* find no larger value. This is a numerical check, not a proof, and the paragraph is
commentary: the theorem in §4 does not depend on it. Hence **1.19281555 is the largest floor8 this certificate chain
(multinomial lower bound for |U−U|, Chernoff bound for |U+U|, symmetric N_h = 0 types) can give at k = 2 with this
mask and base.** Sharper counting of U+U, or using h-blocks with an orientation rule, is outside this statement.
k = 3 at x would give ≈ 1.19359533* (exact Z₃ path-sum computation; asymptotic only, not certified finitely, not claimed).

## 6. Reproduce without downloading (Comments-column summary)
mask ⟨1971,2016,2100,2628,2688,2800⟩∩[0,26972], q = 27022, k = 2, x = 0x1.ffde1bafccc06p-1,
ℓ = 10^18 (2×10^18 digits), N_h = 0, N_r = ⌊ℓ x^{c₂(r)}/Z₂(x)⌋, n₀ = 3362808316958,
C = 26096544438312721237750, |M| = 5869, |M+M| = 31274, c₂ = min-cost residue costs mod q²
(26540 distinct values; generating rule §1), θ > 1.19281555; limit 1.19281555778…*.
