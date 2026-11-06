"""Certificat autonome : 7 points, 21 distances Polybrutus, rang affine 5."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations

K = [
 [0,4,4,12,12,4,4],
 [4,0,4,12,13,9,12],
 [4,4,0,16,9,9,12],
 [12,12,16,0,7,7,12],
 [12,13,9,7,0,13,16],
 [4,9,9,7,13,0,4],
 [4,12,12,12,16,4,0],
]
B = [[K[0][i]+K[0][j]-K[i][j] for j in range(1,7)]
     for i in range(1,7)]  # B = 8 Gram, avec K = 4 d^2

def det(a):
    A = [[F(v) for v in row] for row in a]
    n = len(A); d = F(1)
    for j in range(n):
        pivot = next((i for i in range(j,n) if A[i][j]), None)
        if pivot is None: return F(0)
        if pivot != j:
            A[j],A[pivot] = A[pivot],A[j]; d = -d
        x = A[j][j]; d *= x
        for i in range(j+1,n):
            t = A[i][j]/x
            for k in range(j+1,n): A[i][k] -= t*A[j][k]
            A[i][j] = 0
    return d

def minor(M,ix):
    return det([[M[i][j] for j in ix] for i in ix])

assert Counter(K[i][j] for i,j in combinations(range(7),2)) == {4:6,7:2,9:3,12:6,13:2,16:2}
assert all(K[i][j] == K[j][i] for i in range(7) for j in range(7))
w = [4,-5,-5,4,4,0]
assert all(sum(B[i][j]*w[j] for j in range(6)) == 0 for i in range(6))
I = [0,1,2,3,5]
A = [[B[i][j] for j in I] for i in I]
assert [minor(A,list(range(j))) for j in range(1,6)] == [8,48,1024,3840,14336]
assert [minor(B,list(ix)) for ix in combinations(range(6),5)] == [0,14336,14336,22400,22400,14336]
assert det(B) == 0
print('PASS: mêmes 21 distances, B positive semi-définie, rang affine exactement 5')
