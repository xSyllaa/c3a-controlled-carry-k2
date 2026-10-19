#!/usr/bin/env python3
"""Replay finite-ell k=2 controlled-carry certificate for C_3a.
Reads data/k2_cost_hist.tsv + data/k2_meta.json produced by finite_k2.

Default certification path: N_h=0 (conservative). Residue h=Q/2 is self-inverse,
so the r↔−r orientation injectivity lemma does not cover h-blocks; skeptics 1–3
flagged this. With N_h=0 every lift is a function of the residue and the
distinct-residue injectivity argument applies.

SECONDARY checker (mpmath loggamma at high precision with small pads; not directed rounding).
Primary is cert_primary.py (MPFR directed + Robbins). Run 2026-10-19-c3a-xopt:
x via --xhex (default x* = 0x1.ffde1bafccc06p-1), claim via --claim; exit 0 iff theta_lo > claim.

Credit: Sylla / xSyllaa.
"""
from __future__ import annotations
import argparse, json, sys
from collections import Counter
from pathlib import Path
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
GENS = (1971, 2016, 2100, 2628, 2688, 2800)

def load():
    meta = json.loads((ROOT / "data/k2_meta.json").read_text())
    hist = {}
    for line in (ROOT / "data/k2_cost_hist.tsv").read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        c, m = map(int, line.split())
        hist[c] = m
    return meta, hist

def x_from_hex(h: str) -> "mp.mpf":
    # 0x1.ffffffffffffpN
    assert h.startswith("0x1.") and "p" in h
    body, exp = h[4:].split("p")
    frac = int(body, 16)
    bits = len(body) * 4
    return (mp.mpf(1) + mp.mpf(frac) / mp.power(2, bits)) * mp.power(2, int(exp))

def build_Pplus(B, gens, x):
    reach = bytearray(B + 1)
    reach[0] = 1
    for n in range(B + 1):
        if reach[n]:
            for g in gens:
                if n + g <= B:
                    reach[n + g] = 1
    M = [n for n in range(B + 1) if reach[n]]
    sum_present = bytearray(2 * B + 1)
    for a in M:
        for b in M:
            sum_present[a + b] = 1
    Pp = mp.mpf(0)
    for s in range(2 * B + 1):
        if sum_present[s]:
            Pp += mp.power(x, s)
    return Pp, len(M), int(sum(sum_present))

XHEX = "0x1.ffde1bafccc06p-1"
CLAIM = "1.19281555"

def certify(ell: int = 10**18, dps: int = 120, nh_zero: bool = True) -> dict:
    """Certify finite-ℓ bound.

    nh_zero=True (default, panel fix): omit h=Q/2 blocks entirely (N_h=0).
    nh_zero=False: legacy N_h=2⌊ℓ p_h/2⌋ path (requires orientation argument).
    """
    meta, hist = load()
    B, q, k = meta["B"], meta["q"], 2
    c_h_val = meta["ch"]
    hist_work = dict(hist)
    hist_work[c_h_val] = hist[c_h_val] - 1
    assert hist[0] == 1 and hist_work[c_h_val] % 2 == 0

    mp.mp.dps = dps
    x = x_from_hex(XHEX)
    Z = mp.mpf(0)
    for c, m in hist.items():
        Z += m * mp.power(x, c)
    Pp, nM, nS = build_Pplus(B, GENS, x)

    cnt_mult = Counter()
    used = 0
    C = 0
    ph = mp.power(x, c_h_val) / Z
    if nh_zero:
        nh = 0
    else:
        nh = 2 * int(mp.floor(ell * ph / 2))
        cnt_mult[nh] += 1
        used += nh
        C += nh * c_h_val
    for c, m in hist_work.items():
        if c == 0 or m == 0:
            continue
        n = int(mp.floor(ell * mp.power(x, c) / Z))
        cnt_mult[n] += m
        used += m * n
        C += m * n * c
    n0 = ell - used
    if n0 < 0:
        return {"ok": False, "FAIL": f"n0={n0}", "nh_zero": nh_zero}
    cnt_mult[n0] += 1
    if C % 2:
        return {"ok": False, "FAIL": "odd C", "nh_zero": nh_zero}
    L = C // 2
    Ndig = k * ell

    with mp.workdps(dps + 40):
        logD = mp.loggamma(ell + 1)
        for nval, mult in cnt_mult.items():
            if nval > 0 and mult:
                logD -= mult * mp.loggamma(nval + 1)
        logD_lo = logD - mp.mpf("1e-30") * (1 + sum(cnt_mult.values()))
        logS_hi = Ndig * mp.log(Pp) - C * mp.log(x) + mp.mpf("1e-30")
        logden_hi = mp.log(mp.mpf(4 * B) / mp.mpf(q - 1) * mp.power(mp.mpf(q), Ndig) + 1)
        numer_lo = logD_lo - logS_hi
        if numer_lo <= 0:
            return {"ok": False, "FAIL": "numer_lo<=0", "nh_zero": nh_zero}
        theta_lo = 1 + numer_lo / logden_hi
        # Integer floor8/12 — compare inside high-dps block (leaving workdps
        # can make th8 mpf compare 1e-122 under the decimal literal).
        floor8_int = int(mp.floor(theta_lo * mp.mpf("1e8")))
        floor12_int = int(mp.floor(theta_lo * mp.mpf("1e12")))
        th8_str = f"{floor8_int / 1e8:.8f}"
        th12_str = f"{floor12_int / 1e12:.12f}"
        pass_gt_236 = bool(theta_lo > mp.mpf("1.19230254"))
        pass_gt_page = bool(theta_lo > mp.mpf("1.19102809"))
        pass_gt_claim = bool(theta_lo > mp.mpf(CLAIM))
        pass_gt_retired = bool(theta_lo > mp.mpf("1.19279920"))
        theta_lo_str = mp.nstr(theta_lo, 120, strip_zeros=False)

    return {
        "ok": True,
        "construction": "controlled-carry finite-ell Prop.controlled-block k=2",
        "certification_mode": "N_h=0 conservative (panel fix)" if nh_zero else "legacy N_h>0 (needs h-orientation)",
        "nh_zero": nh_zero,
        "nh": int(nh),
        "c_h": int(c_h_val),
        "credit": "Sylla / xSyllaa",
        "ell": ell,
        "N_digits": Ndig,
        "B": B,
        "q": q,
        "k": 2,
        "gens": list(GENS),
        "|M|": nM,
        "|M+M|": nS,
        "n0": int(n0),
        "C": int(C),
        "L": int(L),
        "sqrtZ2": mp.nstr(mp.sqrt(Z), 30),
        "Pplus": mp.nstr(Pp, 30),
        "theta_lo": theta_lo_str,
        "theta_floor8": th8_str,
        "theta_floor12": th12_str,
        "PASS_gt_1.19230254": pass_gt_236,
        "PASS_gt_1.19102809": pass_gt_page,
        "claim": CLAIM,
        "x_hex": XHEX,
        "PASS_theta_gt_claim": pass_gt_claim,
        "PASS_gt_retired_1.19279920": pass_gt_retired,
        "hypotheses": [
            "0 in M, M=<gens>∩[0,B], B in M",
            "k=2 blocks, base q=27022, residue alphabet Z/q^k Z",
            "c_k min-cost per residue; Z_k = sum x^{c_k(r)}",
            (
                "type counts: N_r=floor(ell p_r) paired r↔−r; N_h=0 (omit self-inverse h=Q/2); N_0=remainder"
                if nh_zero
                else "type counts: N_r=floor(ell p_r) paired; N_h=2*floor(ell p_h/2); N_0=remainder"
            ),
            (
                "injective block encoding via residue-determined lifts (N_h=0 ⇒ distinct-residue lemma applies)"
                if nh_zero
                else "injective oriented block encoding (Prop controlled-block; h-orientation needs prefix rule)"
            ),
            "|U-U| >= multinomial(ell; N_r)",
            "|U+U| <= P_+(x)^{k ell} x^{-C}  (P_+ = sum_{s in M+M} x^s, not pressure root)",
            "denom <= 4 B q^{k ell}/(q-1) + 1 (4maxU+1, V=2U)",
            "GHR-dilated as in K2026b Prop controlled-block",
        ],
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ell", type=int, default=10**18)
    ap.add_argument("--xhex", default="0x1.ffde1bafccc06p-1")
    ap.add_argument("--claim", default="1.19281555")
    ap.add_argument("--dps", type=int, default=120)
    ap.add_argument("--json-out", type=str)
    ap.add_argument(
        "--legacy-nh",
        action="store_true",
        help="Use legacy N_h=2⌊ℓ p_h/2⌋ (needs h-orientation argument). Default is N_h=0.",
    )
    args = ap.parse_args()
    global XHEX, CLAIM
    XHEX, CLAIM = args.xhex, args.claim
    res = certify(args.ell, args.dps, nh_zero=not args.legacy_nh)
    print(json.dumps(res, indent=2))
    if args.json_out:
        Path(args.json_out).write_text(json.dumps(res, indent=2) + "\n")
    ok = (
        res.get("ok")
        and res.get("PASS_gt_1.19230254")
        and res.get("PASS_theta_gt_claim")
    )
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
