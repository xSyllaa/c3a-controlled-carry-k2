#!/usr/bin/env python3
"""PRIMARY certificate: finite controlled-carry construction, block length k (default 3; k=2 reproduces PR #238), N_h=0.

    C_3a  >=  1 + (log|U-U| - log|U+U|) / log(4 maxU + 1)   [GHR lemma with V = 2U]

with the bounds (see THEOREM note):
    |U-U| >= multinomial(ell; N_r)                    (injective type class, N_h = 0)
    |U+U| <= P_+(x)^{k ell} x^{-C}                    (any 0 < x < 1)
    4 maxU + 1 <= 4 B q^{k ell}/(q-1) + 1

All real arithmetic uses MPFR with directed rounding (RNDD/RNDU) through ctypes, ln(n!) via
Robbins (1955) strict bounds  n ln n - n + ln(2 pi n)/2 + 1/(12n+1) < ln n! < ... + 1/(12n).
Integer data (mask, M+M, type counts N_r, n0, C) are exact; N_r = floor(ell x^c / Z_k(x)) is
evaluated from an MPFR enclosure and the run ABORTS if the enclosure straddles an integer.
Requires only python3 (stdlib) + libmpfr.so.6 / libgmp.so.10.

Usage: python3 cert_k3_primary.py --hist ../data/k3_cost_hist.tsv          (k=3 defaults: ell=1e24, x3, claim 1.19359680)
       python3 cert_k3_primary.py --k 2 --hist ../../data/k2_cost_hist.tsv --ell 1000000000000000000 \
           --xhex 0x1.ffde1bafccc06p-1 --claim 1.19281555     (regression: reproduces the k=2 certificate)
Exit code 0 iff  theta_lo > claim  and every structural assertion holds.
Credit: Sylla / xSyllaa.  Based on hostile skeptic #2 (c3a5) MPFR/Robbins checker.
"""
from __future__ import annotations
import argparse, hashlib, json, sys
from decimal import Decimal
from pathlib import Path
from mpfr_ctypes import Mpfr, MPFR_RNDD, MPFR_RNDU, MPFR_RNDN

PREC = 512
B, Q, K = 26972, 27022, 3
QK = Q ** K
GENS = (1971, 2016, 2100, 2628, 2688, 2800)
HIST_SHA256_BY_K = {2: "4d16010f1b41d875c2912f397e76e8d4010e49b635e139ad1f9bcbd741956c63",
                    3: "3186f6e220fd4a36369bf23ce841a9d54e6179cfaff00b5f3189c058f0d85c63"}
HIST_SHA256 = None  # set from --k in main()
PI_LO = "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196"
PI_HI = "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038197"
REFERENCE_BARS = {"PR238_k2": "1.19281555", "PR236": "1.19230254", "page_finite_K2026b": "1.19102809"}
PAGE_LIMIT = "1.19519192"

def M_(v=None): return Mpfr(PREC) if v is None else Mpfr(PREC).set_str(str(v), MPFR_RNDN)

def build_mask_and_S():
    inM = bytearray(B + 1); inM[0] = 1
    for n in range(1, B + 1):
        if any(n >= g and inM[n - g] for g in GENS): inM[n] = 1
    M = [n for n in range(B + 1) if inM[n]]
    S = sorted({a + b for a in M for b in M})
    return M, S

def x_from_hex(h: str) -> tuple[Mpfr, int, int]:
    """exact dyadic x = (2^(4n)+frac) * 2^(e-4n); returns (Mpfr, numerator, log2 denominator)."""
    assert h.startswith("0x1.") and "p" in h, h
    body, e = h[4:].split("p"); e = int(e); n = 4 * len(body)
    num = (1 << n) + int(body, 16); den_log2 = n - e
    x = Mpfr(PREC).set_str(str(num), MPFR_RNDN)
    x.div(x, Mpfr(PREC).set_str(str(1 << den_log2), MPFR_RNDN), MPFR_RNDN)  # exact (power of 2)
    assert num < (1 << den_log2), "need x < 1"
    return x, num, den_log2

_LN2PI = None
def ln2pi():
    global _LN2PI
    if _LN2PI is None:
        two = M_(2)
        lo = Mpfr(PREC).log(Mpfr(PREC).set_str(PI_LO, MPFR_RNDD), MPFR_RNDD); lo.add(lo, Mpfr(PREC).log(two, MPFR_RNDD), MPFR_RNDD)
        hi = Mpfr(PREC).log(Mpfr(PREC).set_str(PI_HI, MPFR_RNDU), MPFR_RNDU); hi.add(hi, Mpfr(PREC).log(two, MPFR_RNDU), MPFR_RNDU)
        _LN2PI = (lo, hi)
    return _LN2PI

def lnfact(n: int, upper: bool) -> Mpfr:
    """Robbins: strict lower (upper=False) / upper (upper=True) bound for ln(n!)."""
    if n <= 1: return Mpfr(PREC).set_si(0)
    r = MPFR_RNDU if upper else MPFR_RNDD
    nd = M_(n)                                   # exact (n < 2^PREC)
    lnn = Mpfr(PREC).log(nd, r)
    t = Mpfr(PREC).mul(Mpfr(PREC).add(nd, M_("0.5"), MPFR_RNDN), lnn, r)  # (n+1/2) ln n
    t.sub(t, nd, r)
    half_l2p = Mpfr(PREC).mul(M_("0.5"), ln2pi()[1 if upper else 0], r)
    t.add(t, half_l2p, r)
    pad = Mpfr(PREC).ui_div(1, M_(12 * n if upper else 12 * n + 1), r)
    t.add(t, pad, r)
    return t

def floor_from_enclosure(lo: Mpfr, hi: Mpfr) -> int:
    a = Decimal(lo.to_str(80, MPFR_RNDD)); b = Decimal(hi.to_str(80, MPFR_RNDU))
    fa, fb = int(a.to_integral_value(rounding="ROUND_FLOOR")), int(b.to_integral_value(rounding="ROUND_FLOOR"))
    if fa != fb: raise SystemExit(f"ABORT: floor enclosure straddles integer [{a},{b}]")
    return fa

def gt(a: Mpfr, b: Mpfr) -> bool:
    return Mpfr(PREC).sub(a, b, MPFR_RNDD).get_d() > 0

def c_h_direct(M):
    """c_k(h), h = q^k/2: carry rule over digits of h, independent of the histogram."""
    INF = 10**9
    kap = {}
    for a in M:
        for b in M:
            d = a - b
            if a + b < kap.get(d, INF): kap[d] = a + b
    h = QK // 2; digs = [(h // Q**i) % Q for i in range(K)]
    best = INF
    from itertools import product
    for gs in product((0, -1), repeat=K):          # g_1..g_K, g_0 = 0
        g = (0,) + gs; tot = 0
        for i in range(K):
            d = digs[i] - g[i] + Q * g[i + 1]
            tot += kap.get(d, INF)
        best = min(best, tot)
    return best

def certify(ell: int, xhex: str, claim: str, hist_path: Path) -> dict:
    raw = hist_path.read_bytes(); sha = hashlib.sha256(raw).hexdigest()
    assert HIST_SHA256 is None or sha == HIST_SHA256, f"histogram sha256 mismatch {sha}"
    hist = [tuple(map(int, l.split())) for l in raw.decode().splitlines() if l.strip() and not l.startswith("#")]
    H = dict(hist)
    assert sum(H.values()) == QK, "histogram must cover all q^k residues"
    assert H[0] == 1
    M, S = build_mask_and_S()
    global CH_EXPECTED
    CH_EXPECTED = c_h_direct(M)
    odd = sorted(c for c, m in hist if m % 2)
    assert len(odd) == 2 and odd[0] == 0, f"exactly 0 and h must be unpaired, got {odd}"
    CH = odd[1]
    assert CH == CH_EXPECTED, f"c(h) from histogram {CH} != direct carry-rule value {CH_EXPECTED}"
    assert len(M) == 5869 and M[0] == 0 and M[-1] == B and len(S) == 31274
    x, xnum, xden2 = x_from_hex(xhex)

    # Z_k(x) enclosure
    Zlo, Zhi = Mpfr(PREC).set_si(0), Mpfr(PREC).set_si(0)
    xc = {}
    for c, m in hist:
        lo = Mpfr(PREC).pow_ui(x, c, MPFR_RNDD); hi = Mpfr(PREC).pow_ui(x, c, MPFR_RNDU)
        xc[c] = (lo, hi)
        Zlo.add(Zlo, Mpfr(PREC).mul_ui(lo, m, MPFR_RNDD), MPFR_RNDD)
        Zhi.add(Zhi, Mpfr(PREC).mul_ui(hi, m, MPFR_RNDU), MPFR_RNDU)
    # P_+(x) upper bound
    Phi = Mpfr(PREC).set_si(0); Plo = Mpfr(PREC).set_si(0)
    for s in S:
        Phi.add(Phi, Mpfr(PREC).pow_ui(x, s, MPFR_RNDU), MPFR_RNDU)
        Plo.add(Plo, Mpfr(PREC).pow_ui(x, s, MPFR_RNDD), MPFR_RNDD)

    ellm = M_(ell)
    used = C = 0; mults = []   # (N_r, multiplicity)
    for c, m in hist:
        if c == 0: continue
        mm = m - 1 if c == CH else m          # N_h = 0: drop h, keep its 2j paired companions
        assert mm % 2 == 0
        lo = Mpfr(PREC).div(Mpfr(PREC).mul(ellm, xc[c][0], MPFR_RNDD), Zhi, MPFR_RNDD)
        hi = Mpfr(PREC).div(Mpfr(PREC).mul(ellm, xc[c][1], MPFR_RNDU), Zlo, MPFR_RNDU)
        n = floor_from_enclosure(lo, hi)
        if n and mm:
            mults.append((n, mm)); used += n * mm; C += n * mm * c
    n0 = ell - used
    assert n0 >= 0 and C % 2 == 0
    L = C // 2

    # log multinomial lower bound
    logD = lnfact(ell, upper=False)
    for n, mm in mults:
        logD.sub(logD, Mpfr(PREC).mul_ui(lnfact(n, True), mm, MPFR_RNDU), MPFR_RNDD)
    logD.sub(logD, lnfact(n0, True), MPFR_RNDD)
    # log|U+U| upper:  N log P+_hi - C log x  (log x < 0: use log x rounded down)
    Nd = K * ell; Ndm = M_(Nd); Cm = M_(C)
    logS = Mpfr(PREC).mul(Ndm, Mpfr(PREC).log(Phi, MPFR_RNDU), MPFR_RNDU)
    logS.sub(logS, Mpfr(PREC).mul(Cm, Mpfr(PREC).log(x, MPFR_RNDD), MPFR_RNDD), MPFR_RNDU)
    # log denominator upper:  log(4B q^N/(q-1) + 1) <= N log q + log(4B/(q-1)) + 1/(4B q^N/(q-1)) ; pad 1e-100
    pre = Mpfr(PREC).div(M_(4 * B), M_(Q - 1), MPFR_RNDU)
    logDen = Mpfr(PREC).mul(Ndm, Mpfr(PREC).log(M_(Q), MPFR_RNDU), MPFR_RNDU)
    logDen.add(logDen, Mpfr(PREC).log(pre, MPFR_RNDU), MPFR_RNDU)
    logDen.add(logDen, Mpfr(PREC).set_str("1e-100", MPFR_RNDU), MPFR_RNDU)
    num = Mpfr(PREC).sub(logD, logS, MPFR_RNDD)
    assert num.get_d() > 0
    theta = Mpfr(PREC).div(num, logDen, MPFR_RNDD); theta.add_ui(theta, 1, MPFR_RNDD)

    # block-length-k asymptote at this x (informational, lower/upper enclosure)
    lq_hi = Mpfr(PREC).log(M_(Q), MPFR_RNDU); lq_lo = Mpfr(PREC).log(M_(Q), MPFR_RNDD)
    a_lo = Mpfr(PREC).sub(Mpfr(PREC).div_ui(Mpfr(PREC).log(Zlo, MPFR_RNDD), K, MPFR_RNDD), Mpfr(PREC).log(Phi, MPFR_RNDU), MPFR_RNDD)
    a_hi = Mpfr(PREC).sub(Mpfr(PREC).div_ui(Mpfr(PREC).log(Zhi, MPFR_RNDU), K, MPFR_RNDU), Mpfr(PREC).log(Plo, MPFR_RNDD), MPFR_RNDU)
    a_lo.div(a_lo, lq_hi, MPFR_RNDD); a_lo.add_ui(a_lo, 1, MPFR_RNDD)
    a_hi.div(a_hi, lq_lo, MPFR_RNDU); a_hi.add_ui(a_hi, 1, MPFR_RNDU)

    th_s = theta.to_str(60, MPFR_RNDD)
    a, b = th_s.split(".")
    claim_m = Mpfr(PREC).set_str(claim, MPFR_RNDU)
    out = {
        "checker": "cert_k3_primary.py (MPFR directed rounding + Robbins), credit Sylla / xSyllaa",
        "construction": f"finite controlled-carry, K2026b Prop.3 instantiated at k={K}, N_h=0, V=2U",
        "hist_sha256": sha, "B": B, "q": Q, "k": K, "gens": list(GENS),
        "|M|": len(M), "|M+M|": len(S), "distinct_costs": len(hist), "c_h": CH,
        "x_hex": xhex, "x_rational": f"{xnum}/2^{xden2}",
        "ell": str(ell), "N_digits": str(Nd), "N_h": 0, "n0": str(n0), "C": str(C), "L": str(L),
        "nonzero_type_classes": sum(mm for _, mm in mults),
        "Zk_lo": Zlo.to_str(30, MPFR_RNDD), 
        "Pplus_hi": Phi.to_str(30, MPFR_RNDU),
        "logD_lo": logD.to_str(40, MPFR_RNDD), "logS_hi": logS.to_str(40, MPFR_RNDU), "logDen_hi": logDen.to_str(40, MPFR_RNDU),
        "theta_lo": th_s, "theta_floor8": f"{a}.{b[:8]}", "theta_floor12": f"{a}.{b[:12]}",
        "k_asymptote_at_x": [a_lo.to_str(20, MPFR_RNDD), a_hi.to_str(20, MPFR_RNDU)],
        "claim": claim, "PASS_theta_gt_claim": gt(theta, claim_m),
        **{f"PASS_gt_{k}_{v}": gt(theta, Mpfr(PREC).set_str(v, MPFR_RNDU)) for k, v in REFERENCE_BARS.items()},
        "PASS_below_page_limit_1.19519192": gt(Mpfr(PREC).set_str(PAGE_LIMIT, MPFR_RNDD), theta),
    }
    out["PASS_all"] = all(v for k, v in out.items() if k.startswith("PASS_"))
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=3, choices=(2, 3))
    ap.add_argument("--ell", type=int, default=10**24)
    ap.add_argument("--xhex", default="0x1.ffde3d0d5025dp-1")
    ap.add_argument("--claim", default="1.19359680")
    ap.add_argument("--hist", required=True)
    ap.add_argument("--json-out")
    a = ap.parse_args()
    global K, QK, HIST_SHA256
    K = a.k; QK = Q ** K; HIST_SHA256 = HIST_SHA256_BY_K[K]
    res = certify(a.ell, a.xhex, a.claim, Path(a.hist))
    s = json.dumps(res, indent=2); print(s)
    if a.json_out: Path(a.json_out).write_text(s + "\n")
    print(("PASS" if res["PASS_all"] else "FAIL") + f": theta_lo = {res['theta_floor12']}... > {a.claim} : {res['PASS_theta_gt_claim']}")
    sys.exit(0 if res["PASS_all"] else 1)

if __name__ == "__main__":
    main()
