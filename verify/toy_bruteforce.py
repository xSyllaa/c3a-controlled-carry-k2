#!/usr/bin/env python3
"""Toy brute force of the EXACT N_h=0 construction (even q, k=2, residue-determined min-cost lifts,
r <-> -r paired types, L = C/2, V = 2U semantics via 4maxU+1).  For tiny q (4, 6) it enumerates
U, U-U, U+U exactly and asserts every lemma used by the big certificate:
  injectivity (distinct u-v over the whole type class), u,v in U (weight = L), |U-U| >= multinomial,
  |U+U| <= P_+^{2 ell} x^{-C}, and theta_bound <= theta_exact.
stdlib only.  Origin: hostile skeptic #4 (c3a4) toy.py; sympy dependency removed.
Credit: Sylla / xSyllaa."""
import itertools, math
def multiset_permutations(letters):
    return set(itertools.permutations(letters))
def run(q,M,ell,x,pairs):
    Q=q*q; B=max(M)
    ka={}
    for a in M:
        for b in M:
            d=a-b
            if d not in ka or a+b<ka[d][0]: ka[d]=(a+b,a,b)
    # residue-determined lift per residue mod Q: min cost e=d0+d1*q
    lift={}
    for r in range(Q):
        best=None
        for d1 in ka:
            for d0 in ((r-d1*q)%Q,(r-d1*q)%Q-Q):
                if d0 in ka:
                    c=ka[d0][0]+ka[d1][0]
                    if best is None or c<best[0]: best=(c,d0,d1)
        assert best, r
        lift[r]=best
    # enforce antisymmetry lift(-r) = -lift(r) for chosen representative
    for r in range(1,Q):
        nr=(-r)%Q
        if r<nr: lift[nr]=(lift[r][0],-lift[r][1],-lift[r][2])
    assert all(lift[r][0]==lift[(-r)%Q][0] for r in range(Q))
    h=Q//2
    Z=sum(x**lift[r][0] for r in range(Q))
    Nr={}
    for r,n in pairs:
        assert r!=h and r!=0
        Nr[r]=n; Nr[(-r)%Q]=n
    n0=ell-sum(Nr.values()); C=sum(Nr[r]*lift[r][0] for r in Nr); L=C//2
    assert C%2==0
    letters=[0]*n0+[r for r in Nr for _ in range(Nr[r])]
    def real(r):  # (a-digits, b-digits) for the block: use ka for r rep, swap for -r
        c,d0,d1=lift[r]
        a0,b0=(ka[d0][1],ka[d0][2]) if d0>=0 else (ka[-d0][2],ka[-d0][1])
        a1,b1=(ka[d1][1],ka[d1][2]) if d1>=0 else (ka[-d1][2],ka[-d1][1])
        return (a0,a1),(b0,b1)
    # check swap consistency: realization for -r is swap of r
    U=set()
    for digs in itertools.product(M,repeat=2*ell):
        if sum(digs)<=L: U.add(sum(d*q**i for i,d in enumerate(digs)))
    ds=set(); cnt=0
    for seq in multiset_permutations(letters):
        u=v=0
        for i,r in enumerate(seq):
            (a0,a1),(b0,b1)=real(r)
            assert a0-b0+(a1-b1)*q==lift[r][1]+lift[r][2]*q
            u+=(a0+a1*q)*Q**i; v+=(b0+b1*q)*Q**i
        assert u in U and v in U, "u/v not in U"
        ds.add(u-v); cnt+=1
    UmU={a-b for a in U for b in U}; UpU={a+b for a in U for b in U}
    P=sum(x**s for s in {a+b for a in M for b in M})
    Sbound=P**(2*ell)*x**(-C)
    th=1+math.log(len(UmU)/len(UpU))/math.log(2*max(U)+1)
    thb=1+(math.log(cnt)-math.log(Sbound))/math.log(4*max(U)+1)
    print(f"q={q} M={M} ell={ell} type_size={cnt} distinct_d={len(ds)} |U|={len(U)} |U-U|={len(UmU)} |U+U|={len(UpU)} Sbound={Sbound:.1f} theta_exact={th:.5f} theta_bound={thb:.5f}")
    assert len(ds)==cnt and len(UmU)>=cnt and len(UpU)<=Sbound+1e-9 and thb<=th+1e-12
run(4,[0,1,3],4,0.8,[(1,1),(5,1)])
run(4,[0,1,3],4,0.8,[(3,1),(6,1)])
run(4,[0,1,3],5,0.8,[(1,1),(4,1)])
run(6,[0,1,2,4],3,0.8,[(7,1)])
run(6,[0,1,3,4],4,0.8,[(5,1),(13,1)])
print("TOY OK")
