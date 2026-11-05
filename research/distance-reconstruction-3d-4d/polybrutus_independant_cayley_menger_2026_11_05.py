#!/usr/bin/env python3
"""2e enumeration INDEPENDANTE (Cayley-Menger, pas de trilateration
par Gram du premier enumerateur) des distances non etiquetees en R3.
Mode laboratoire sans ecriture, Python standard, entiers exacts.
"""
from collections import Counter
from itertools import combinations, permutations, product
import time
import laboratoire_spectral as geometry

K=(4,7,9,12,13,16)
MAX=Counter({4:6,7:2,9:3,12:6,13:2,16:2})
EDGES=tuple(combinations(range(7),2))
REFERENCE=tuple(tuple(0 if i==j else int(4*geometry.d2e(i,j)) for j in range(7)) for i in range(7))
assert Counter(REFERENCE[i][j] for i,j in EDGES)==MAX
REFERENCE_SIGNATURES={
    tuple(REFERENCE[p[i]][p[j]] for i,j in EDGES)
    for p in permutations(range(7))
}
assert len(REFERENCE_SIGNATURES)==2520

def det_int(M):
    """Elimination fraction-free sur des entiers, sans flottants."""
    a=[list(map(int,row)) for row in M]
    n=len(a)
    sign=1;den=1
    for h in range(n-1):
        v=next((i for i in range(h,n) if a[i][h]!=0),None)
        if v is None:return 0
        if v!=h:a[h],a[v]=a[v],a[h];sign=-sign
        pivot=a[h][h]
        for i in range(h+1,n):
            for j in range(h+1,n):
                num=a[i][j]*pivot-a[i][h]*a[h][j]
                assert num%den==0
                a[i][j]=num//den
        for i in range(h+1,n):a[i][h]=0
        den=pivot
    return sign*a[n-1][n-1]

def cm(D,indices):
    n=len(indices)
    M=[[0]+[1]*n]
    for i in indices:
        M.append([1]+[D[i][j] for j in indices])
    return det_int(M)

def allowed(used):
    return all(used[k]<=MAX[k] for k in K)

def root_matrix(edges,base):
    D=[[0]*base for _ in range(base)]
    for (i,j),w in zip(combinations(range(base),2),edges):
        D[i][j]=D[j][i]=w
    return D

def root_key(D):
    n=len(D)
    return min(tuple(D[p[i]][p[j]] for i,j in combinations(range(n),2))
               for p in permutations(range(n)))

def root_valid(D):
    n=len(D)
    for size in range(3,n+1):
        sign=(-1)**size
        if any(sign*cm(D,group)<=0 for group in combinations(range(n),size)):
            return False
    return True

def extend_once(D,w):
    n=len(D)
    E=[row+[w[i]] for i,row in enumerate(D)]
    E.append(list(w)+[0])
    return E

def pair_allowed(D,wa,wb,weight,base):
    # Le determinant CM de base+1 points doit s'annuler sur chaque
    # sous-ensemble de (base-1) racines et les 2 candidats.
    M=[row[:]+[wa[i],wb[i]] for i,row in enumerate(D)]
    M.append(list(wa)+[0,weight])
    M.append(list(wb)+[weight,0])
    a,b=base,base+1
    return all(cm(M,subset+(a,b))==0
               for subset in combinations(range(base),base-1))

def final_matrix(D,chosen,links,base):
    T=[row[:]+[w[i] for w in chosen] for i,row in enumerate(D)]
    for i,wi in enumerate(chosen):
        row=list(wi)
        for j in range(len(chosen)):
            row.append(0 if i==j else links[min(i,j),max(i,j)])
        T.append(row)
    return T

def enumerate_rank(rank,stop_when_different=True):
    base=rank+1
    seen=set();counts=Counter();different=[]
    begin=time.monotonic()
    rootedges=tuple(combinations(range(base),2))
    for values in product(K,repeat=len(rootedges)):
        usedroot=Counter(values)
        if not allowed(usedroot):continue
        D=root_matrix(values,base)
        key=root_key(D)
        if key in seen:continue
        seen.add(key)
        if not root_valid(D):continue
        counts['valid_roots']+=1
        # Toutes les distances possibles entre chaque nouveau point et racines.
        options=[]
        for w in product(K,repeat=base):
            cc=Counter(w)
            if not allowed(usedroot+cc):continue
            if cm(extend_once(D,w),tuple(range(base+1)))!=0:continue
            options.append((w,cc))
        if len(options)<7-base:continue
        counts['root_with_enough_extensions']+=1
        counts['possible_single_extensions']+=len(options)
        used=usedroot.copy()
        chosen=[]
        links={}

        def dfs(start):
            if len(chosen)==7-base:
                T=final_matrix(D,chosen,links,base)
                sig=tuple(T[i][j] for i,j in EDGES)
                assert Counter(sig)==MAX
                # Tous les determinants correspondant a rank+2 points doivent
                # etre nuls si la dimension affine est au plus rank.
                assert all(cm(T,subset)==0 for subset in combinations(range(7),base+1))
                counts['complete']+=1
                if sig not in REFERENCE_SIGNATURES:
                    counts['noncongruent']+=1
                    if len(different)<2:different.append(T)
                    return stop_when_different
                counts['congruent']+=1
                return False
            remaining=7-base-len(chosen)
            for ix in range(start,len(options)-remaining+1):
                w,cc=options[ix]
                if not allowed(used+cc):continue
                possibilities=[]
                for prev in chosen:
                    ks=[k for k in K if allowed(used+cc+Counter({k:1}))
                        and pair_allowed(D,prev,w,k,base)]
                    if not ks:break
                    possibilities.append(ks)
                else:
                    for assignment in product(*possibilities):
                        inc=cc.copy()
                        inc.update(assignment)
                        if not allowed(used+inc):continue
                        n=len(chosen)
                        for j,k in enumerate(assignment):links[j,n]=k
                        chosen.append(w); used.update(inc)
                        if dfs(ix+1):return True
                        chosen.pop();used.subtract(inc)
                        for j in range(n):del links[j,n]
            return False
        if dfs(0):break
    counts['canonical_roots_examined']=len(seen)
    print('RANK',rank,dict(counts),'seconds',round(time.monotonic()-begin,2),flush=True)
    if different:print('NONCONGRUENT WITNESS',different[0])
    return counts,different

def main():
    for rank in (3,2,1):
        counts,different=enumerate_rank(rank)
        assert not different,'Uniqueness disproved!'
        if rank==3:assert counts['complete']>0
        else:assert counts['complete']==0
    print('PASS - INDEPENDENT CAYLEY-MENGER ENUMERATION in R3')
    print('NO R3 NONCONGRUENT CONFIGURATION IN EXACT DISTANCE MULTISET')
    print('METHOD: full root edge assignments; CM determinants; no import of first search algorithm')
    print('PRODUCTION_WRITES=0')

if __name__=='__main__':
    main()