#!/usr/bin/env python3
"""Verification mathematique independante des matrices finales de la recherche
sans contrainte de sphere, avec SymPy et Gram exactes. LAB ONLY.
"""
from collections import Counter
from itertools import combinations,permutations
import sympy as S
import polybrutus_unicite_sans_sphere_2026_11_05 as scan

EDGES=tuple(combinations(range(7),2))
REF=S.Matrix(scan.REF)
ALL_SIGNATURES=set(tuple(scan.REF[p[i]][p[j]] for i,j in EDGES)
                   for p in permutations(range(7)))
assert len(ALL_SIGNATURES)==2520
checks=[]
orig=scan.signature

def checked_signature(base,edges,chosen,weights):
    seq,M=orig(base,edges,chosen,weights)
    assert len(seq)==21
    assert Counter(seq)==scan.COUNTS
    D=S.Matrix(M)
    assert D==D.T
    assert all(D[i,i]==0 for i in range(7))
    assert all(D[i,j]>0 for i,j in EDGES)
    # Reconstruire la matrice de Gram ancree en P0, sans hypothese de sphere.
    G=S.Matrix([[(D[0,i]+D[0,j]-D[i,j])/S.Integer(2)
                 for j in range(1,7)] for i in range(1,7)])
    # Une matrice de Gram semidefinie-positive de rang <=3 est embeddable en 3D.
    assert G.rank()==base-1,(base,G.rank())
    for order in range(1,4):
        for sub in combinations(range(6),order):
            minor=G.extract(sub,sub).det(method='bareiss')
            assert minor>=0,(base,order,sub,minor)
    # Independamment de l'ancrage, une configuration non etiquetee
    # equivalente doit fournir le rang D=4 (configuration cospherique).
    assert D.rank()==4
    assert seq in ALL_SIGNATURES
    checks.append((base,tuple(seq)))
    return seq,M

scan.signature=checked_signature

def main():
    for base in (4,3,2):
        tallies,ex=scan.search_rank(base,stop_on_counterexample=False)
        assert not ex
        if base==4:assert tallies['valid_configurations']==27
        else:assert tallies['valid_configurations']==0
    assert len(checks)==27 and all(base==4 for base,_ in checks)
    print('SYMPY PASS: 27 matrices finales independamment controlees')
    print('SYMPY PASS: 21 distances, Gram PSD exact, rang affine 3, rang D=4')
    print('SYMPY PASS: chaque configuration correspond par permutation a la source')
    print('SYMPY PASS: aucune solution plane ou alignee')
    print('CONCLUSION: les memes 21 distances non etiquetees forcent la meme sphere de rayon 1 a isometrie pres')
    print('LIMITE: sous realisabilite dans R3 et distances exactes uniquement.')
    print('PRODUCTION_WRITES=0')
if __name__=='__main__':
    main()