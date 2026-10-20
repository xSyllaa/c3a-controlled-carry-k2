"""Independent skeptic replay for the k=3 controlled-carry finite bound (fresh code path).
- mask, kappa, M+M rebuilt with numpy (no shared code with builders / cert_k3_primary.py)
- histogram: total, symmetry parity, c(h) by brute force over all lifts, exact low-cost prefix by
  direct enumeration of digit triples, sampled residues by brute force
- theta: mpmath interval arithmetic (mpmath.iv), Robbins bounds, floor straddle check
usage: skeptic_k3.py HIST ELL XHEX CLAIM
"""
import sys, hashlib, json, random
import numpy as np
from mpmath import iv, mpf, mp
G=(1971,2016,2100,2628,2688,2800); B=26972; q=27022; K=3; QK=q**K; h=QK//2
hist_path, ELL, XHEX, CLAIM = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]

# --- mask via boolean dynamic programming
inM=np.zeros(B+1,bool); inM[0]=True
for n in range(1,B+1):
    for g in G:
        if n>=g and inM[n-g]: inM[n]=True; break
M=np.flatnonzero(inM); assert len(M)==5869
# kappa(d) = min_{a in M, a-d in M} 2a-d, d in [-B,B]; compute by looping over a
INF=10**9
kap=np.full(2*B+1,INF,np.int64)          # index d+B
for a in M:
    b=M                                   # d=a-b
    idx=a-b+B; val=a+b
    np.minimum.at(kap,idx,val)
Dmask=kap<INF; assert Dmask.sum()==52193
assert kap[B]==0
assert all(kap[B+d]==kap[B-d] for d in range(B+1))
S=np.unique(np.add.outer(M,M)); assert len(S)==31274
def kp(d):  # kappa or INF
    return int(kap[d+B]) if -B<=d<=B else INF
def c_brute(r):
    """min over ALL (d0,d1,d2) in D^3 with d0+d1 q+d2 q^2 = r mod q^3, |d_i|<=B<q.
    d0 = r mod q or that minus q; then recurse. Enumerates all candidates explicitly."""
    best=INF
    r%=QK
    for d0 in (r%q, r%q-q):
        k0=kp(d0)
        if k0>=INF: continue
        r1=((r-d0)//q)%(q*q)
        for d1 in (r1%q, r1%q-q):
            k1=kp(d1)
            if k1>=INF: continue
            r2=((r1-d1)//q)%q
            for d2 in (r2, r2-q):
                k2=kp(d2)
                if k2<INF: best=min(best,k0+k1+k2)
    return best

raw=open(hist_path,'rb').read(); sha=hashlib.sha256(raw).hexdigest()
H={}
for l in raw.decode().splitlines():
    if l.startswith('#') or not l.strip(): continue
    c,m=map(int,l.split()); assert c not in H; H[c]=m
assert sum(H.values())==QK
odd=sorted(c for c,m in H.items() if m%2)
ch=c_brute(h); assert odd==[0,ch], (odd,ch)
assert H[0]==1 and ch==24367
# exact low-cost prefix: residues with c<=T enumerated from triples with kappa-sum<=T
T=7000
small=[d for d in range(-B,B+1) if kap[d+B]<=T]
best={}
sk=[(d,int(kap[d+B])) for d in small]
for d0,k0 in sk:
    for d1,k1 in sk:
        if k0+k1>T: continue
        for d2,k2 in sk:
            s=k0+k1+k2
            if s>T: continue
            r=(d0+d1*q+d2*q*q)%QK
            if best.get(r,INF)>s: best[r]=s
pref={}
for r,s in best.items(): pref[s]=pref.get(s,0)+1
for c in range(T+1):
    assert pref.get(c,0)==H.get(c,0), ("prefix mismatch",c,pref.get(c,0),H.get(c,0))
nprefix=sum(pref.values())
# sampled residues: brute-force costs must lie in histogram support; compare empirical mean
random.seed(20261020); samp=[c_brute(random.randrange(QK)) for _ in range(200000)]
assert all(s in H for s in samp)
emp=sum(samp)/len(samp)
hm=sum(c*m for c,m in H.items())/QK
hv=sum(c*c*m for c,m in H.items())/QK-hm*hm
z=(emp-hm)/(hv/len(samp))**.5

# --- theta in interval arithmetic
iv.prec=600; mp.prec=600
body,e=XHEX[4:].split('p'); e=int(e); nb=4*len(body)
num=(1<<nb)+int(body,16); x=iv.mpf(num)/iv.mpf(2)**(nb-int(e)) if False else iv.mpf(num)*iv.mpf(2)**(e-nb)
assert x.a==x.b and x<1
lx=iv.log(x)
def xp(c): return iv.exp(lx*c)
Z=iv.mpf(0)
pw={}
for c,m in H.items():
    pw[c]=xp(c); Z+=pw[c]*m
P=iv.mpf(0)
for s in S: P+=xp(int(s))
ell=iv.mpf(ELL)
used=0; C=0; terms=[]
for c,m in H.items():
    if c==0: continue
    mm=m-1 if c==ch else m
    v=ell*pw[c]/Z
    fa,fb=int(mp.floor(v.a)),int(mp.floor(v.b))
    assert fa==fb, ("floor straddle",c)
    if fa and mm: terms.append((fa,mm)); used+=fa*mm; C+=fa*mm*c
n0=ELL-used; assert n0>=0 and C%2==0
L2PI=iv.log(2*iv.pi)
def lnf(n):   # Robbins: ln n! in [s+1/(12n+1), s+1/(12n)]
    if n<2: return iv.mpf(0)
    nn=iv.mpf(n); s=(nn+iv.mpf(0.5))*iv.log(nn)-nn+L2PI/2
    return iv.mpf([ (s+1/(12*nn+1)).a, (s+1/(12*nn)).b ])
logM=lnf(ELL)-lnf(n0)
for n,mm in terms: logM-=lnf(n)*mm
N=K*ELL
logS=iv.log(P)*N-lx*C
logDen=iv.log(iv.mpf(4*B)/(q-1))+iv.log(iv.mpf(q))*N+iv.mpf('1e-100')
th=1+(logM-logS)/logDen
thlo=th.a
claim=iv.mpf(CLAIM)
out=dict(sha=sha,c_h=ch,prefix_T=T,prefix_residues=nprefix,sample_n=len(samp),sample_mean=emp,hist_mean=hm,z=z,
 n0=str(n0),C=str(C),Z_lo=mp.nstr(Z.a,25),P_hi=mp.nstr(P.b,25),theta_lo=mp.nstr(thlo,30),theta_width=mp.nstr(th.b-th.a,5),
 PASS=bool(thlo>claim.b), PASS_gt_1_19281555=bool(thlo>iv.mpf('1.19281555').b))
print(json.dumps(out,indent=1))
