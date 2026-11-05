#!/usr/bin/env python3
"""Recherche exhaustive (rangs affines 1,2,3) des 21 distances non etiquetees.
Distance carree D_ij = k_ij/4, k entiers; placements dans R^3, SANS spheres.
Programme laboratoire: calculs entiers exacts, aucune ecriture en production.
"""
from collections import Counter
from itertools import combinations,permutations,product
import time
import laboratoire_spectral as lab

VALUES=(4,7,9,12,13,16)
COUNTS=Counter({4:6,7:2,9:3,12:6,13:2,16:2})
PAIRS=tuple(combinations(range(7),2))
REF=tuple(tuple(0 if i==j else int(4*lab.d2e(i,j)) for j in range(7)) for i in range(7))
assert Counter(REF[i][j] for i,j in PAIRS)==COUNTS
FINGERPRINTS=set()
for perm in permutations(range(7)):
    FINGERPRINTS.add(tuple(REF[perm[i]][perm[j]] for i,j in PAIRS))
print('Empreintes isometriques precompute:',len(FINGERPRINTS),flush=True)

def can_sub(have,need):
    return all(need[k]<=have[k] for k in VALUES)

def gram_data(anchor,base):
    # anchor: base edge weights in order combinations(base,2)
    m=[[0]*base for _ in range(base)]
    for (i,j),k in zip(combinations(range(base),2),anchor):
        m[i][j]=m[j][i]=k
    H=[[m[0][i]+m[0][j]-m[i][j] for j in range(1,base)] for i in range(1,base)]
    r=base-1
    if r==1:
        a=H[0][0]
        return (a,((1,),)) if a>0 else None
    if r==2:
        a,u=H[0];_,b=H[1]
        det=a*b-u*u
        return (det,((b,-u),(-u,a))) if a>0 and det>0 else None
    a,u,v=H[0];_,b,w=H[1];_,_,c=H[2]
    d=a*b*c+2*u*v*w-a*w*w-b*v*v-c*u*u
    if a<=0 or a*b-u*u<=0 or d<=0:return None
    return d,((b*c-w*w,v*w-u*c,u*w-b*v),
              (v*w-u*c,a*c-v*v,u*v-a*w),
              (u*w-b*v,u*v-a*w,a*b-u*u))

def dot(A,x,y):
    return sum(x[i]*A[i][j]*y[j] for i in range(len(x)) for j in range(len(y)))

def canonical_anchor(arr,base):
    mat=[[0]*base for _ in range(base)]
    for (i,j),k in zip(combinations(range(base),2),arr):
        mat[i][j]=mat[j][i]=k
    return min(tuple(mat[p[i]][p[j]] for i,j in combinations(range(base),2))
               for p in permutations(range(base)))

def candidate_list(base,anchors,left,det,adj):
    n=base-1
    pairs=list(combinations(range(base),2))
    m=[[0]*base for _ in range(base)]
    for (i,j),k in zip(pairs,anchors):
        m[i][j]=m[j][i]=k
    candidates=[]
    for dist in product(VALUES,repeat=base):
        req=Counter(dist)
        if not can_sub(left,req):continue
        q0=dist[0]
        v=tuple(m[0][i]+q0-dist[i] for i in range(1,base))
        if dot(adj,v,v)==2*q0*det:
            candidates.append((dist,v,req))
    return candidates

def pairweight(A,B,det,adj):
    a,va,_=A; b,vb,_=B
    num=(a[0]+b[0])*det-dot(adj,va,vb)
    if num%det:return None
    k=num//det
    return k if k in COUNTS else None

def signature(base,edges,chosen,weights):
    M=[[0]*7 for _ in range(7)]
    for (i,j),k in zip(combinations(range(base),2),edges):
        M[i][j]=M[j][i]=k
    for q,ind in enumerate(chosen):
        ix=base+q
        for i,k in enumerate(ind[0]):
            M[i][ix]=M[ix][i]=k
    for (i,j),k in weights.items():
        M[base+i][base+j]=M[base+j][base+i]=k
    return tuple(M[i][j] for i,j in PAIRS),M

def search_rank(base,stop_on_counterexample=True):
    start=time.monotonic()
    canonical_seen=set()
    tally=Counter()
    counterexamples=[]
    cases=0
    root_perms=tuple(combinations(range(base),2))
    for anchor in product(VALUES,repeat=len(root_perms)):
        req=Counter(anchor)
        if not can_sub(COUNTS,req):continue
        key=canonical_anchor(anchor,base)
        if key in canonical_seen:continue
        canonical_seen.add(key)
        gd=gram_data(anchor,base)
        if gd is None:continue
        det,adj=gd
        tally['admissible_anchor']+=1
        left=COUNTS-req
        cand=candidate_list(base,anchor,left,det,adj)
        if len(cand)<7-base:continue
        tally['candidate_total']+=len(cand)
        pairmap={}
        for i in range(len(cand)):
            for j in range(i+1,len(cand)):
                k=pairweight(cand[i],cand[j],det,adj)
                if k is not None:pairmap[i,j]=k
        tally['compatible_candidate_pairs']+=len(pairmap)
        selected=[]
        used=Counter()
        pairweights={}
        def rec(pos):
            if len(selected)==7-base:
                tally['valid_configurations']+=1
                chosen=[cand[i] for i in selected]
                seq,_=signature(base,anchor,chosen,pairweights)
                if seq not in FINGERPRINTS:
                    tally['noncongruent']+=1
                    if len(counterexamples)<3: counterexamples.append(dict(anchor=anchor,
                        new_vertex_anchor_distances=[list(c[0]) for c in chosen],
                        interpoint=list(pairweights.values()),signature=seq))
                    return stop_on_counterexample
                tally['congruent']+=1
                return False
            need=(7-base)-len(selected)
            if len(cand)-pos<need:return False
            for i in range(pos,len(cand)-need+1):
                candidate=cand[i]
                changes=Counter(candidate[2])
                temp=[]
                okay=True
                for j,sel in enumerate(selected):
                    ij=(sel,i)
                    k=pairmap.get(ij)
                    if k is None:okay=False;break
                    temp.append((len(selected),j,k))
                    changes[k]+=1
                if not okay or not can_sub(left-used,changes):continue
                used.update(changes); selected.append(i)
                for a,b,k in temp:
                    pairweights[b,a]=k
                stop=rec(i+1)
                for a,b,_ in temp:del pairweights[b,a]
                selected.pop()
                used.subtract(changes)
                if stop:return True
            return False
        stopped=rec(0)
        if stopped:break
    tally['canonical_anchors_total']=len(canonical_seen)
    tally['elapsed_seconds']=round(time.monotonic()-start,3)
    print('RANG',base-1,dict(tally),flush=True)
    if counterexamples:
        print('CONTRE_EXEMPLE EXACT',repr(counterexamples[0]),flush=True)
    return tally,counterexamples

def main():
    for base in (4,3,2):
        tally,examples=search_rank(base,stop_on_counterexample=True)
        if examples:
            print('UNICITE GLOBALE REFUTEE (rang affine '+str(base-1)+')')
            print('LAB_ONLY_NO_PRODUCTION_WRITES')
            return
    print('PASS: aucun contre-exemple avec rang affine 1,2,3 en R3')
    print('Unicite non etiquetee en R3 verifiee PAR ENUMERATION ENTIERE sous distances exactes.')
    print('LIMITES: rayon non fixe; distances non perturbees; nouveaute scientifique non etablie.')
    print('LAB_ONLY_NO_PRODUCTION_WRITES')

if __name__=='__main__':main()