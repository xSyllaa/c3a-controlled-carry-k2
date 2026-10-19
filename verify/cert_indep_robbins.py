# INDEPENDENT checker (hostile skeptic #4 lineage): k=2, N_h=0, Robbins-Stirling bounds, mpmath 80 dps.
# Reads its OWN histogram/sumset built by hist_indep.c (run hist_indep in a scratch dir first).
# Usage: python3 cert_indep_robbins.py ELL XHEX HIST.tsv SUMSET.txt
# Credit: Sylla / xSyllaa.
import sys, json
from mpmath import mp, mpf, log, pi, floor, nstr
mp.dps = 80
ell = int(sys.argv[1]) if len(sys.argv)>1 else 10**18
xh = sys.argv[2] if len(sys.argv)>2 else '0x1.ffde1bafccc06p-1'
B, q = 26972, 27022; Q2 = q*q; ch = 24367
x = mpf(float.fromhex(xh))
hist = {}
for line in open(sys.argv[3] if len(sys.argv)>3 else 'my_hist.tsv'):
    if line.startswith('#'): continue
    c, m = map(int, line.split()); hist[c] = m
assert sum(hist.values()) == Q2 and hist[0] == 1
Z2 = sum(m * x**c for c, m in hist.items())
# exclude r=0 and h (self-inverse); N_h = 0
counts = []  # (cost, multiplicity, N_r)
for c, m in hist.items():
    if c == 0: continue
    mult = m - (1 if c == ch else 0)
    if mult == 0: continue
    Nr = int(floor(ell * x**c / Z2))
    counts.append((c, mult, Nr))
n0 = ell - sum(m*N for c, m, N in counts)
C = sum(m*N*c for c, m, N in counts)
assert n0 >= 0 and C % 2 == 0
L = C // 2
def lf_lo(n):  # lower bound log n!
    if n == 0: return mpf(0)
    n = mpf(n); return n*log(n) - n + log(2*pi*n)/2 + 1/(12*n+1)
def lf_hi(n):
    if n == 0: return mpf(0)
    n = mpf(n); return n*log(n) - n + log(2*pi*n)/2 + 1/(12*n)
logD = lf_lo(ell) - lf_hi(n0) - sum(m*lf_hi(N) for c, m, N in counts)
S = [int(l) for l in open(sys.argv[4] if len(sys.argv)>4 else 'sumset.txt')]
assert len(S) == 31274
Pp = sum(x**s for s in S)
N = 2*ell
logS = N*log(Pp) - C*log(x)
logDen = log(mpf(4*B)*mpf(q)**N/(q-1) + 1)
theta = 1 + (logD - logS)/logDen
# lemma side: V=2U gives |V-V|<=maxV+1<2maxV+1 automatically; also check D<=den sanity
out = dict(ell=ell, N_digits=N, nh=0, n0=n0, C=C, L=L, Z2=nstr(Z2,25), sqrtZ2=nstr(mp.sqrt(Z2),20),
           Pplus=nstr(Pp,20), logD=nstr(logD,25), logS=nstr(logS,25), logDen=nstr(logDen,25),
           theta=nstr(theta,30), floor8=nstr(floor(theta*10**8)/10**8,10),
           gt_1_19281555=bool(theta>mpf('1.19281555')), gt_1_19279920=bool(theta>mpf('1.19279920')),
           gt_236=bool(theta>mpf('1.19230254')),
           asym=nstr(1+(log(Z2)-2*log(Pp))/(2*log(q)), 15))
print(json.dumps(out, indent=1))
