#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Polybrutus–Higgs: exact exploratory geometry, no physical prediction.

The input 3D/4D distance data are taken from the dated Polybrutus research
certificates.  Arithmetic is Fraction-exact; Python standard library only.

The standard electroweak Higgs potential is discussed analytically but NOT
claimed to generate any seven-point geometry or any boson-mass prediction.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from typing import Sequence

N = 7
EXPECTED = Counter({4: 6, 7: 2, 9: 3, 12: 6, 13: 2, 16: 2})

# Coordinates (sqrt(weight[j]) * entry[j]); exactly the reference source.
POLY3 = (
    (0, 0, 1),
    (F(1, 2), 0, F(1, 2)),
    (F(-1, 4), F(3, 4), F(1, 2)),
    (F(-1, 4), F(-3, 4), F(1, 2)),
    (0, 0, -1),
    (F(1, 2), 0, F(-1, 2)),
    (F(-1, 2), 0, F(-1, 2)),
)
W3 = (3, 1, 1)

Q4_A = (
    (0, 0, 0, 0),
    (1, 0, 0, 0),
    (F(1, 2), F(1, 2), 0, 0),
    (F(1, 2), F(-1, 6), F(2, 3), 0),
    (F(3, 8), F(11, 24), F(13, 24), F(1, 8)),
    (F(-1, 8), F(-1, 24), F(7, 24), F(-1, 8)),
    (1, F(1, 3), F(2, 3), 0),
)
W4_A = (1, 3, 6, 30)

# B and C: exact distance certificates recovered locally, K_ij=4d_ij^2.
K4_B = (
    (0, 9, 16, 12, 12, 4, 16),
    (9, 0, 12, 9, 4, 9, 13),
    (16, 12, 0, 12, 7, 13, 4),
    (12, 9, 12, 0, 4, 12, 4),
    (12, 4, 7, 4, 0, 7, 4),
    (4, 9, 13, 12, 7, 0, 12),
    (16, 13, 4, 4, 4, 12, 0),
)
K4_C = (
    (0, 12, 13, 4, 9, 7, 9),
    (12, 0, 4, 16, 9, 4, 12),
    (13, 4, 0, 12, 13, 12, 4),
    (4, 16, 12, 0, 7, 12, 4),
    (9, 9, 13, 7, 0, 4, 12),
    (7, 4, 12, 12, 4, 0, 16),
    (9, 12, 4, 4, 12, 16, 0),
)

# Six cube-face normals + reference on S^3, theta=60 deg.
CUBE7 = (
    (1, 0, 0, 0),
    (F(1, 2), F(1, 2), 0, 0),
    (F(1, 2), F(-1, 2), 0, 0),
    (F(1, 2), 0, F(1, 2), 0),
    (F(1, 2), 0, F(-1, 2), 0),
    (F(1, 2), 0, 0, F(1, 2)),
    (F(1, 2), 0, 0, F(-1, 2)),
)
W_CUBE = (1, 3, 3, 3)


def matrix_from_coordinates(coords, weights):
    """Return K=4 d^2 with rational entries."""
    return tuple(tuple(
        sum((4 * F(w) * (F(a) - F(b)) ** 2
             for w, a, b in zip(weights, coords[i], coords[j])), F(0))
        for j in range(len(coords))) for i in range(len(coords)))


def histogram(K):
    return Counter(F(K[i][j]) for i, j in combinations(range(len(K)), 2))


def inverse(M):
    n = len(M)
    X = [[F(M[i][j]) for j in range(n)]
         + [F(int(i == j)) for j in range(n)] for i in range(n)]
    for k in range(n):
        pivot = next((i for i in range(k, n) if X[i][k]), None)
        if pivot is None:
            raise ValueError('singular matrix')
        X[k], X[pivot] = X[pivot], X[k]
        d = X[k][k]
        X[k] = [q / d for q in X[k]]
        for i in range(n):
            if i != k and X[i][k]:
                t = X[i][k]
                X[i] = [a - t * b for a, b in zip(X[i], X[k])]
    return [row[n:] for row in X]


def rank(M):
    X = [[F(t) for t in row] for row in M]
    if not X:
        return 0
    m, n = len(X), len(X[0])
    r = 0
    for j in range(n):
        i = next((i for i in range(r, m) if X[i][j]), None)
        if i is None:
            continue
        X[r], X[i] = X[i], X[r]
        d = X[r][j]
        X[r] = [a / d for a in X[r]]
        for i in range(m):
            if i != r and X[i][j]:
                t = X[i][j]
                X[i] = [a - t * b for a, b in zip(X[i], X[r])]
        r += 1
        if r == m:
            break
    return r



def determinant(M):
    X = [[F(t) for t in row] for row in M]
    n = len(X)
    d = F(1)
    for j in range(n):
        r = next((k for k in range(j, n) if X[k][j]), None)
        if r is None:
            return F(0)
        if r != j:
            X[j], X[r] = X[r], X[j]
            d = -d
        pivot = X[j][j]
        d *= pivot
        for k in range(j+1, n):
            t = X[k][j] / pivot
            X[k] = [a - t * b for a,b in zip(X[k],X[j])]
    return d


def certify_euclidean(K, dim):
    """Exact PSD proof: positive definite dim-anchors + vanishing Schur.

    Establishes existence of a Euclidean realization of K in R^dim.
    """
    G = anchored_gram(K)
    A = [r[:dim] for r in G[:dim]]
    minors = tuple(determinant([r[:j] for r in A[:j]])
                   for j in range(1,dim+1))
    assert all(t > 0 for t in minors), minors
    Ai = inverse(A)
    for i in range(dim, N - 1):
        for j in range(dim, N - 1):
            predicted = sum((G[i][a] * Ai[a][b] * G[b][j]
                             for a in range(dim) for b in range(dim)), F(0))
            assert predicted == G[i][j], ('nonzero Schur remainder',i,j)
    assert rank(G) == dim
    return minors


def anchored_gram(K):
    return [[F(K[0][i] + K[0][j] - K[i][j], 8)
             for j in range(1, N)] for i in range(1, N)]


def centered_gram(K):
    """-1/8 H K H, where K_ij=4*d_ij^2, H=I-J/7."""
    row = [sum(map(F, r), F(0)) for r in K]
    total = sum(row, F(0))
    return [[-F(K[i][j], 8) + (row[i] + row[j]) / 56 - total / 392
             for j in range(N)] for i in range(N)]


def sphere_certificate(K, dim):
    """Exact affine-sphere test using dim anchor vectors and Gram identities.

    Returns sphere center norm^2 and residuals; for dim=4, the first five
    points must be affinely independent (all four sample matrices satisfy).
    """
    G = anchored_gram(K)
    A = [r[:dim] for r in G[:dim]]
    Ai = inverse(A)
    b = [G[i][i] / 2 for i in range(dim)]
    coeff = [sum((Ai[i][j] * b[j] for j in range(dim)), F(0))
             for i in range(dim)]
    center_sq = sum((coeff[i] * b[i] for i in range(dim)), F(0))
    residual = [sum((2 * G[i][j] * coeff[j] for j in range(dim)), F(0))
                - G[i][i] for i in range(N - 1)]
    return {
        'cospherical': all(t == 0 for t in residual),
        'radius_sq': center_sq,
        'residual': tuple(residual),
    }


def pair_energy(K, inverse_power=1):
    """E_2=sum (d_ij^2)^(-inverse_power), for nonzero distances."""
    return sum(((F(4, K[i][j])) ** inverse_power
                for i, j in combinations(range(N), 2)), F(0))


def three_body_energy(K):
    """E_3=sum_{i<j<k} (d_ij^2)(d_ik^2)(d_jk^2); a descriptor only."""
    return sum((F(K[i][j] * K[i][k] * K[j][k], 64)
                for i, j, k in combinations(range(N), 3)), F(0))


def equilateral_counts(K):
    return Counter(F(K[i][j], 4)
                   for i, j, k in combinations(range(N), 3)
                   if K[i][j] == K[i][k] == K[j][k])


def centered_moment_3(K):
    G = centered_gram(K)
    return sum((G[i][j] * G[j][k] * G[k][i]
                for i in range(N) for j in range(N) for k in range(N)), F(0))


def tangent_force_sq(K, R2=F(1)):
    """Squared projected gradient of sum(1/d_ij^2) on S^(d-1).

    Uses only the unit-vector Gram u_i.u_j=1-d_ij^2/(2R^2);
    not an experimentally motivated interaction or a gauge invariant.
    """
    U = [[F(1) - F(K[i][j], 8 * R2) for j in range(N)]
         for i in range(N)]
    total = F(0)
    for i in range(N):
        # Cartesian force: -2 sum_{j!=i}(u_i-u_j)/d_ij^4.
        c = [F(0) for _ in range(N)]
        for j in range(N):
            if i == j:
                continue
            a = F(-2) / (F(K[i][j], 4) ** 2)
            c[i] += a
            c[j] -= a
        norm2 = sum((c[j] * c[k] * U[j][k]
                     for j in range(N) for k in range(N)), F(0))
        radial = sum((c[j] * U[i][j] for j in range(N)), F(0))
        total += norm2 - radial ** 2
    assert total >= 0
    return total


def run():
    K3 = matrix_from_coordinates(POLY3, W3)
    A = matrix_from_coordinates(Q4_A, W4_A)
    cubes = matrix_from_coordinates(CUBE7, W_CUBE)
    inputs = [('Polybrutus_3D', K3, 3), ('4D_A', A, 4),
              ('4D_B', K4_B, 4), ('4D_C', K4_C, 4)]
    details = {}
    for name, K, dim in inputs:
        assert len(K) == N and all(len(row) == N for row in K)
        assert all(K[i][i] == 0 and K[i][j] == K[j][i]
                   for i in range(N) for j in range(N))
        assert histogram(K) == EXPECTED, (name, histogram(K))
        assert sum(K[i][j] for i, j in combinations(range(N), 2)) == 195
        assert rank(centered_gram(K)) == dim
        minors = certify_euclidean(K, dim)
        sph = sphere_certificate(K, dim)
        if name != '4D_B':
            assert sph['cospherical'] and sph['radius_sq'] == 1, (name, sph)
        else:
            assert not sph['cospherical']
        details[name] = {
            'rank': dim,
            'positive_anchor_minors': minors,
            'spherical': sph['cospherical'],
            'sphere_radius_sq_anchor': sph['radius_sq'],
            'E2': pair_energy(K),
            'E3': three_body_energy(K),
            'moment3': centered_moment_3(K),
            'equi': equilateral_counts(K),
            'projected_force_norm_sq': tangent_force_sq(K) if sph['cospherical'] else None,
        }

    assert all(v['E2'] == F(6329, 546) for v in details.values())
    assert [details[k]['E3'] for k in ('Polybrutus_3D', '4D_A', '4D_B', '4D_C')] == [
        F(24435, 64), F(23835, 64), F(26817, 64), F(24273, 64)]

    assert rank(centered_gram(cubes)) == 4
    assert histogram(cubes) == Counter({4: 6, 6: 12, 12: 3})
    assert sum(CUBE7[i][0] for i in range(7)) == 4
    assert pair_energy(cubes) == 15

    # Exact potential V(lambda,v,h)=lambda/4 * (2 v h+h^2)^2.
    # At h=0: V=V'=0 and V''=2 lambda v^2.
    # Tangent unit-vector directions have zero V for all finite angular changes.
    return details, {
        'cube_barycenter_e0': F(4, 7),
        'cube_E2': pair_energy(cubes),
        'cube_E3': three_body_energy(cubes),
        'higgs_radial_mass_squared_over_lambda_v_squared': 2,
        'higgs_angular_modes_with_zero_potential_curvature': 3,
    }


def format_fraction(q):
    if q is None:
        return 'NON APPLICABLE'
    return str(q)


if __name__ == '__main__':
    cases, control = run()
    print('POLYBRUTUS–HIGGS | EXPLORATION EXACTE | PASS')
    for name, d in cases.items():
        print(f'{name}: rank={d["rank"]}, sphere={d["spherical"]}, '
              f'E2={d["E2"]}, E3={d["E3"]}, '
              f'trGram3={d["moment3"]}, '
              f'projected_F2={format_fraction(d["projected_force_norm_sq"])}')
    print('Cube+1:', control)
    print('PHYSIQUE: aucune masse prédite; jauge non résolue; pas de preuve Higgs.')
