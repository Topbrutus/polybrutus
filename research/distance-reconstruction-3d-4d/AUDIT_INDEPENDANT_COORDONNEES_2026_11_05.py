#!/usr/bin/env python3
"""Self-contained exact-coordinate audit for external mathematical reviewers.
No import from Polybrutus libraries; only SymPy. No output files, no network.
"""
from collections import Counter
from itertools import combinations, permutations
import sympy as S

R = S.Rational
a,b,c = S.sqrt(3),S.sqrt(6),S.sqrt(30)
P=(
 (0,0,1),
 (a/2,0,R(1,2)),
 (-a/4,R(3,4),R(1,2)),
 (-a/4,-R(3,4),R(1,2)),
 (0,0,-1),
 (a/2,0,-R(1,2)),
 (-a/2,0,-R(1,2))
)
Q=(
 (0,0,0,0),
 (1,0,0,0),
 (R(1,2),a/2,0,0),
 (R(1,2),-a/6,2*b/3,0),
 (R(3,8),11*a/24,13*b/24,c/8),
 (-R(1,8),-a/24,7*b/24,-c/8),
 (1,a/3,2*b/3,0)
)
edges=tuple(combinations(range(7),2))
expected=Counter({R(1):6,R(7,4):2,R(9,4):3,R(3):6,R(13,4):2,R(4):2})

def squared(x,y):
    return S.simplify(sum(((u-v)**2 for u,v in zip(x,y)),S.Integer(0)))

def analyze(points):
    D=tuple(tuple(S.Integer(0) if i==j else squared(points[i],points[j])
                  for j in range(7)) for i in range(7))
    hist=Counter(D[i][j] for i,j in edges)
    v=[S.Matrix(points[i])-S.Matrix(points[0]) for i in range(1,7)]
    gram=S.Matrix([[S.simplify(x.dot(y)) for y in v] for x in v])
    return D,hist,gram

def main():
    dp,hp,gp=analyze(P)
    dq,hq,gq=analyze(Q)
    assert hp==expected and hq==expected,(hp,hq,expected)
    assert gp.rank()==3 and gq.rank()==4
    assert S.factor(gq[:4,:4].det())==R(15,16)
    assert S.simplify(S.sqrt(gq[:4,:4].det())/S.factorial(4)) == S.sqrt(15)/96
    # A complete comparison of every relabeling rules out a hidden congruence
    # (equivalent labeled distance matrices).
    equivalents=0
    for perm in permutations(range(7)):
        if all(dp[i][j]==dq[perm[i]][perm[j]] for i,j in edges):
            equivalents+=1
    assert equivalents==0
    # Positive control: order changes cannot change the multiset.
    assert analyze(list(reversed(Q)))[1]==expected
    # Negative control: a non-isometric perturbation fails.
    shifted=list(Q)
    shifted[4]=tuple(shifted[4][j]+(R(1,17) if j==3 else 0) for j in range(4))
    assert analyze(shifted)[1]!=expected
    print('PASS: canonical P and witness Q have the same exact 21-distance multiset')
    print('PASS: all 5040 permutations checked, matching labeled matrices =',equivalents)
    print('PASS: affine Gram ranks P=3, Q=4')
    print('PASS: positive 4-simplex Gram determinant = 15/16; volume = sqrt(15)/96')
    print('PASS: positive permutation and negative perturbation controls')
    print('NOTE: this audits the explicit witness only; it does not independently enumerate R3.')
    print('PRODUCTION_WRITES=0')

if __name__=='__main__':
    main()