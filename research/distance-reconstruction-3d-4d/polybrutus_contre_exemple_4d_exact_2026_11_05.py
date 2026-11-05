#!/usr/bin/env python3
"""Contre-exemple 4D exact aux 21 distances NON ETIQUETEES de Polybrutus.
Ce script n'utilise PAS le chercheur de configurations: il verifie des
coordonnees explicites avec SymPy. Aucun fichier de production modifie.
"""
import sympy as S
from itertools import combinations
from collections import Counter

sq3,sq6,sq30=S.sqrt(3),S.sqrt(6),S.sqrt(30)
Q=S.Rational
POINTS=(
  (0,0,0,0),
  (1,0,0,0),
  (Q(1,2),sq3/2,0,0),
  (Q(1,2),-sq3/6,2*sq6/3,0),
  (Q(3,8),11*sq3/24,13*sq6/24,sq30/8),
  (-Q(1,8),-sq3/24,7*sq6/24,-sq30/8),
  (1,sq3/3,2*sq6/3,0),
)
TARGET=Counter({Q(1):6,Q(7,4):2,Q(9,4):3,Q(3):6,Q(13,4):2,Q(4):2})
EXPECTED_K=(
 (0,4,4,12,12,4,16),
 (4,0,4,12,13,9,12),
 (4,4,0,16,9,9,12),
 (12,12,16,0,7,7,4),
 (12,13,9,7,0,13,4),
 (4,9,9,7,13,0,12),
 (16,12,12,4,4,12,0)
)
def dot(a,b):
    return sum((x*y for x,y in zip(a,b)),S.Integer(0))

def d2(a,b):
    return S.simplify(sum(((x-y)**2 for x,y in zip(a,b)),S.Integer(0)))

def main():
    P=[S.Matrix(v) for v in POINTS]
    assert len(P)==7 and all(len(p)==4 for p in P)
    actual=Counter()
    for i,j in combinations(range(7),2):
        v=d2(P[i],P[j])
        actual[v]+=1
        assert 4*v==EXPECTED_K[i][j],(i,j,v,EXPECTED_K[i][j])
    assert actual==TARGET,(actual,TARGET)
    from laboratoire_spectral import EXACT_POINTS,d2e
    assert Counter(S.Rational(d2e(i,j).numerator,d2e(i,j).denominator) for i,j in combinations(range(7),2))==TARGET
    # Grams with physical distances (not integers k=4 d2)
    differences=[P[i]-P[0] for i in range(1,7)]
    G=S.Matrix([[dot(a,b) for b in differences] for a in differences])
    assert G.rank()==4
    assert G[:4,:4].det()==S.Rational(15,16)
    for order in range(1,7):
        for group in combinations(range(6),order):
            minor=G.extract(group,group).det(method='bareiss')
            assert minor>=0,(order,group,minor)
    assert all(G.extract(g,g).det()==0 for g in combinations(range(6),5))
    # Compare with the original 3D gram matrix, rank = 3.
    old=[S.Matrix([S.sqrt(3)*a,b,c]) for a,b,c in EXACT_POINTS]
    orig=[old[i]-old[0] for i in range(1,7)]
    g0=S.Matrix([[dot(a,b) for b in orig] for a in orig])
    assert g0.rank()==3
    assert G.rank()!=g0.rank()
    assert all(d2(old[i],old[j])==d2e(i,j) for i,j in combinations(range(7),2))
    print('PASS - 7 points 4D explicites, toutes les 21 distances exactes')
    print('FAMILLES:', ', '.join(f'{v}: {actual[v]}' for v in sorted(actual)))
    print('MATRICE K = 4*d² (7x7):')
    for row in EXPECTED_K: print(' ', ' '.join(f'{v:2d}' for v in row))
    print('GRAM rank (contre-exemple):',G.rank())
    print('GRAM rank (Polybrutus canonique):',g0.rank())
    print('GRAM mineur principal 4x4:',G[:4,:4].det())
    print('GRAM PSD: tous les mineurs principaux non-negatifs')
    print('CONCLUSION: deux configurations NON CONGRUENTES en R4, memes distances non etiquetees')
    print('LIMITES: nouveaute dans la litterature non etablie, aucune consequence physique')
    print('PRODUCTION_WRITES=0')

if __name__=='__main__':
    main()