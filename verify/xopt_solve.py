#!/usr/bin/env python3
"""Solve the k=2 x-stationarity condition  g(x) = (1/2) E_p[c] - E_{P+}[s] = 0
(p_r ∝ x^{c(r)} over the cost histogram, P+ weights ∝ x^s over M+M), which maximizes the
k=2 asymptote  theta_inf(x) = 1 + (log Z2(x)/2 - log P+(x)) / log q.
Then compare candidate dyadic x values.  Scouting/optimisation only: the certificate is valid
for ANY x in (0,1); x only needs to be good, not exactly optimal.
Credit: Sylla / xSyllaa."""
import sys
from pathlib import Path
import mpmath as mp
ROOT = Path(__file__).resolve().parents[1]
mp.mp.dps = 60
q = 27022
hist = [tuple(map(int, l.split())) for l in (ROOT/"data/k2_cost_hist.tsv").read_text().splitlines() if l.strip() and not l.startswith('#')]
S = [int(l) for l in (ROOT/"data/sumset.txt").read_text().split()]
def parts(x):
    Z = mp.fsum(m*x**c for c, m in hist); Zc = mp.fsum(m*c*x**c for c, m in hist)
    P = mp.fsum(x**s for s in S); Ps = mp.fsum(s*x**s for s in S)
    return Z, Zc, P, Ps
def th(x):
    Z, _, P, _ = parts(x); return 1 + (mp.log(Z)/2 - mp.log(P))/mp.log(q)
def g(x):
    Z, Zc, P, Ps = parts(x); return Zc/Z/2 - Ps/P
lo, hi = mp.mpf('0.9995'), mp.mpf('0.99999')
assert g(lo) * g(hi) < 0
for _ in range(120):
    mid = (lo+hi)/2
    if (g(mid) > 0) == (g(lo) > 0): lo = mid
    else: hi = mid
xs = (lo+hi)/2
print("x_star (60 dps)      =", mp.nstr(xs, 40))
print("theta_inf(x_star)    =", mp.nstr(th(xs), 25))
def hexval(h):
    body, e = h[4:].split('p'); return (1 + mp.mpf(int(body, 16))/mp.mpf(2)**(4*len(body)))*mp.mpf(2)**int(e)
for h in sys.argv[1:] or ["0x1.ffde827adc0fep-1", "0x1.ffde1bafccc06p-1"]:
    xv = hexval(h)
    print(f"{h:26s} x={mp.nstr(xv,20)}  g={mp.nstr(g(xv),5)}  theta_inf={mp.nstr(th(xv),20)}  loss_vs_opt={mp.nstr(th(xs)-th(xv),5)}")
