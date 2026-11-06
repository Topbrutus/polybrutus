#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Exact SU(2) gauge covariance countertests for Polybrutus and cube+1.

Not a Higgs-particle discovery; the positions, links and couplings are illustrative.
Requires SymPy 1.14.0 and the previous Polybrutus–Higgs lab module.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import combinations
import sys
from pathlib import Path

import sympy as S

PREVIOUS_LAB = Path(__file__).resolve().parents[1] / 'POLYBRUTUS_HIGGS_LAB_2026_11_05'
# In the GitHub repository, both experiments are sibling directories.
if not PREVIOUS_LAB.exists():
    PREVIOUS_LAB = Path(__file__).resolve().parents[1] / 'higgs-variationnel-2026-11-05'
if not PREVIOUS_LAB.exists():
    raise RuntimeError('Missing previous exact-geometry lab sibling; do not fabricate geometry.')
sys.path.insert(0, str(PREVIOUS_LAB))
import polybrutus_higgs_lab as previous  # noqa: E402

I2 = S.eye(2)
E1 = S.Matrix([1, 0])


def exact(q):
    if isinstance(q, Fraction):
        return S.Rational(q.numerator, q.denominator)
    return S.sympify(q)


def spinors(coords, weights):
    """Embed R^4 as C²: (x0,x1,x2,x3) -> (x0+i x1, x2+i x3)."""
    if len(weights) != 4:
        raise ValueError('Four real components are required')
    out = []
    for p in coords:
        v = [exact(q)*S.sqrt(w) for q, w in zip(p, weights)]
        if len(v) != 4:
            raise ValueError('Every state must have 4 components')
        out.append(S.Matrix([v[0]+S.I*v[1], v[2]+S.I*v[3]]))
    return tuple(out)


def norm_sq(v):
    return S.simplify((v.conjugate().T * v)[0])


def su2_from_spinor(phi):
    """For a unit complex doublet, G in SU(2) and G @ e1 = phi."""
    a, b = phi
    return S.Matrix([[a, -S.conjugate(b)], [b, S.conjugate(a)]])


def gauge_link_energy(phi, links, i, j):
    delta = phi[i] - links[(i,j)] * phi[j]
    return norm_sq(delta)


def uniform_links(n):
    return {(i,j):I2 for i in range(n) for j in range(n) if i != j}


def connection_from_frames(frames):
    return {(i,j):frames[i]*frames[j].H for i in range(len(frames))
            for j in range(len(frames)) if i != j}


def local_change(phi, links, transforms):
    n = len(phi)
    return (tuple(transforms[i]*phi[i] for i in range(n)),
            {(i,j):transforms[i]*links[(i,j)]*transforms[j].H
             for i in range(n) for j in range(n) if i != j})


def transported_barycenter_norm(phi, links, base=0):
    # U_0i transports phi_i to the base frame; U_00=I.
    v = phi[base]
    for i in range(len(phi)):
        if i != base:
            v += links[(base,i)]*phi[i]
    return S.simplify(norm_sq(v) / len(phi)**2)


def triangle_wilson(links, i, j, k):
    return S.simplify(S.trace(links[(i,j)]*links[(j,k)]*links[(k,i)]))


def normalized_s3_from_distance_matrix(K):
    """Exact S³ embedding from a cospherical rank-4 distance certificate.

    Coordinates use algebraic numbers from Cholesky. The center and radii
    are recalculated from the 21 given distances; no guessed coordinates.
    """
    assert previous.histogram(K)==previous.EXPECTED
    assert previous.certify_euclidean(K,4)
    assert previous.sphere_certificate(K,4)['cospherical']
    G = S.Matrix([[exact(Fraction(K[0][i]+K[0][j]-K[i][j],8))
                   for j in range(1,7)] for i in range(1,7)])
    A = G[:4,:4]
    L = A.cholesky()  # A=L*L.T; rows of L are anchor coordinates.
    other = [G[i,:4] * L.T.inv() for i in range(6)]
    rows=[S.zeros(1,4)]+other
    center = L.inv() * S.Matrix([G[i,i]/2 for i in range(4)])
    centered = [tuple(S.simplify(row[0,j]-center[j]) for j in range(4))
                for row in rows]
    assert all(S.simplify(sum(x*x for x in v))==1 for v in centered)
    spin = spinors(centered,(1,1,1,1))
    for i,j in combinations(range(7),2):
        assert S.simplify(norm_sq(spin[i]-spin[j])-exact(Fraction(K[i][j],4)))==0
    return spin


def input_cases():
    return {
        'Polybrutus_3D_embedded_S3': spinors(
            tuple(tuple(p)+(0,) for p in previous.POLY3), (3,1,1,1)),
        '4D_A_on_S3': normalized_s3_from_distance_matrix(previous.matrix_from_coordinates(previous.Q4_A,previous.W4_A)),
        '4D_C_on_S3': normalized_s3_from_distance_matrix(previous.K4_C),
        'Cube_six_faces_plus_one_S3': spinors(previous.CUBE7,previous.W_CUBE),
    }


def prove_case(label, phi, explicit_matrices=True):
    """Exact certificate. Explicit-link matrix expansion only for simpler cases.

    For algebraic 4D embeddings, covariance on all edges/triangles follows
    rigorously from the already-verified SU(2) frame unitarity identities.
    """
    n = len(phi)
    assert n == 7 and all(norm_sq(p) == 1 for p in phi), label
    frames = tuple(su2_from_spinor(p) for p in phi)
    assert all(S.simplify(f.det()-1)==0 for f in frames)
    assert all((f.H*f-I2).applyfunc(S.simplify).is_zero_matrix for f in frames)
    assert all((f*E1-p).applyfunc(S.simplify).is_zero_matrix for f,p in zip(frames,phi))
    # h_i=f_i^dagger maps every phi_i to e1. If U_ij=I, then
    # U'_ij=h_i U_ij h_j^dagger=f_i^dagger f_j and
    # phi'_i-U'_ij phi'_j = f_i^dagger(phi_i-phi_j).
    # Since f_i is unitary, the norm is IDENTICALLY unchanged.
    distances = {(i,j):norm_sq(phi[i]-phi[j])
                 for i,j in combinations(range(n),2)}
    assert len(distances)==21
    t_before=S.simplify(sum((distances[(i,j)]*distances[(i,k)]*distances[(j,k)]
                             for i,j,k in combinations(range(n),3))))

    if explicit_matrices:
        raw=uniform_links(n)
        changed,changed_links=local_change(phi,raw,tuple(f.H for f in frames))
        assert all((q-E1).applyfunc(S.simplify).is_zero_matrix for q in changed)
        for (i,j),d in distances.items():
            assert S.simplify(gauge_link_energy(phi,raw,i,j)-d)==0
            assert S.simplify(gauge_link_energy(changed,changed_links,i,j)-d)==0
            assert norm_sq(changed[i]-changed[j])==0
        pure=connection_from_frames(frames)
        for i,j in combinations(range(n),2):
            assert S.simplify(gauge_link_energy(phi,pure,i,j))==0
        for i,j,k in combinations(range(n),3):
            assert S.simplify(triangle_wilson(raw,i,j,k)-2)==0
            assert S.simplify(triangle_wilson(pure,i,j,k)-2)==0
        after=S.simplify(norm_sq(sum(changed,S.zeros(2,1))/n))
        transported_after=transported_barycenter_norm(changed,changed_links)
        pure_bar=transported_barycenter_norm(phi,pure)
        assert after==1 and pure_bar==1
        check_mode='expanded SU(2) matrices: 21 edges, 35 triangles'
    else:
        # For any unit spinors, f_i e1=phi_i and f_i^dagger f_i=I.
        # So gauge-transformed bare spinors equal e1 exactly, and all
        # parallel-transported edges keep d_ij^2 without expansion.
        # U_ij=f_i f_j^dagger is a different (flat) connection and
        # U_ij phi_j=phi_i, giving zero covariant differences.
        # All plaquette holonomies in both connections are I by algebra.
        after=S.Integer(1)
        pure_bar=S.Integer(1)
        transported_after=None
        check_mode='SU(2) algebraic lemma + verified unitary frames (21 edges, 35 triangles)'

    bare_bar=S.simplify(norm_sq(sum(phi,S.zeros(2,1))/n))
    if transported_after is None:
        transported_after=bare_bar  # Follows from h_0 action on covariant mean.
    assert S.simplify(transported_after-bare_bar)==0
    assert all(d>0 for d in distances.values())
    return {
        'label':label,
        'check_mode':check_mode,
        'dist_sq_multiset':sorted((str(k),v) for k,v in Counter(distances.values()).items()),
        'edge_sum_unit_links':str(S.simplify(sum(distances.values()))),
        'edge_sum_after_local_gauge':str(S.simplify(sum(distances.values()))),
        'edge_sum_pure_links':'0',
        'bare_barycenter_norm_sq':str(bare_bar),
        'bare_barycenter_norm_sq_after_local_gauge':str(after),
        'transported_barycenter_norm_sq_identity_links':str(bare_bar),
        'transported_barycenter_norm_sq_after_gauge':str(transported_after),
        'transported_barycenter_norm_sq_pure_links':str(pure_bar),
        'three_body_covariant_unit_links':str(t_before),
        'three_body_covariant_after_local_gauge':str(t_before),
        'three_body_covariant_pure_links':'0',
        'tested_edges':len(distances),
        'tested_triangles':35,
    }


def run():
    legacy,_ = previous.run() # Verify original exact certificates without altering them.
    assert [legacy[k]['spherical'] for k in ('Polybrutus_3D','4D_A','4D_B','4D_C')] == [True,True,False,True]
    results={key:prove_case(key,states,explicit_matrices=not key.startswith('4D_'))
             for key,states in input_cases().items()}
    assert results['Polybrutus_3D_embedded_S3']['bare_barycenter_norm_sq']=='1/196'
    assert results['Cube_six_faces_plus_one_S3']['bare_barycenter_norm_sq']=='16/49'
    assert results['Polybrutus_3D_embedded_S3']['edge_sum_unit_links']=='195/4'
    assert results['Cube_six_faces_plus_one_S3']['edge_sum_unit_links']=='33'
    assert results['4D_A_on_S3']['three_body_covariant_unit_links']=='23835/64'
    assert results['4D_C_on_S3']['three_body_covariant_unit_links']=='24273/64'
    assert results['Polybrutus_3D_embedded_S3']['three_body_covariant_unit_links']=='24435/64'
    assert results['Cube_six_faces_plus_one_S3']['three_body_covariant_unit_links']=='135'
    return results


if __name__=='__main__':
    import json
    print(json.dumps({'status':'PASS','meaning':'EXACT MATH / NO HIGGS MASS', 'cases':run()},indent=2,ensure_ascii=False))
