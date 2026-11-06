#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Seven-site electroweak gauge-covariant STATIC toy model, exact where stated.

No Higgs-mass prediction, no fitted parameter and no continuum 3+1D claim.
Four Polybrutus geometries remain distinguishable mathematical inputs, not
physical evidence. Required: sympy>=1.12 and two earlier sibling labs.
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from fractions import Fraction
from itertools import combinations
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
GEOMETRY_DIR = ROOT / 'higgs-variationnel-2026-11-05'
GAUGE_DIR = ROOT / 'higgs-jauge-2026-11-05'
if not GEOMETRY_DIR.exists():
    GEOMETRY_DIR = ROOT / 'POLYBRUTUS_HIGGS_LAB_2026_11_05'
if not GAUGE_DIR.exists():
    GAUGE_DIR = ROOT / 'POLYBRUTUS_HIGGS_JAUGE_2026_11_05'
if not GEOMETRY_DIR.is_dir() or not GAUGE_DIR.is_dir():
    raise RuntimeError('Missing audited source laboratories; do not synthesize inputs.')
sys.path.insert(0, str(GEOMETRY_DIR))
sys.path.insert(0, str(GAUGE_DIR))
import polybrutus_higgs_lab as geom  # noqa: E402
import higgs_gauge_links as gauge  # noqa: E402

# Parameters fixed BEFORE inspecting any Polybrutus energy; no physical units.
N = 7
A = sp.Integer(1)            # 1D nearest-neighbor arc-length convention
V = sp.sqrt(2)               # Higgs vacuum norm V / sqrt(2) = 1
LAMBDA = sp.Integer(1)
KAPPA = sp.Integer(1)
G_SU2 = sp.Integer(1)        # benchmark couplings, not phenomenology
G_U1 = sp.Integer(1)
Y_HIGGS = sp.Rational(1, 2)

# A cyclic 1D lattice embedded as a regular heptagon in a spatial plane.
# These are space positions; four-component field values live in INTERNAL C^2.
RADIUS = A/(2*sp.sin(sp.pi/N))
EDGES = tuple((i, (i+1) % N) for i in range(N))
I2 = sp.eye(2)
# Predeclared nontrivial flat GLOBAL holonomy on the ring: no local plaquettes.
# SU(2) rotation and a U(1) unit phase in the Higgs representation.
R_SU2 = sp.Matrix([[sp.Rational(3, 5), -sp.Rational(4, 5)],
                   [sp.Rational(4, 5),  sp.Rational(3, 5)]])
Q_U1 = sp.Rational(3, 5) + sp.I*sp.Rational(4, 5)
TWIST = Q_U1 * R_SU2


def check_fixed_protocol():
    assert N == 7 and len(EDGES) == N and A == 1 and KAPPA == 1
    assert LAMBDA == 1 and sp.simplify(V**2) == 2
    assert Y_HIGGS == sp.Rational(1, 2)
    assert (R_SU2.T*R_SU2 - I2) == sp.zeros(2)
    assert R_SU2.det() == 1
    assert sp.expand(Q_U1*sp.conjugate(Q_U1)) == 1
    assert (TWIST.H*TWIST-I2).applyfunc(sp.simplify) == sp.zeros(2)
    assert sp.simplify(2*RADIUS*sp.sin(sp.pi/N)-A) == 0
    return True


def positions():
    return tuple((RADIUS*sp.cos(2*sp.pi*i/N),
                  RADIUS*sp.sin(2*sp.pi*i/N), sp.Integer(0), sp.Integer(0))
                 for i in range(N))


def unit_links(twisted=False):
    return {(i,j): (TWIST if twisted and i==N-1 else I2)
            for i,j in EDGES}


def norm2(vec):
    return sp.simplify((vec.H*vec)[0])


def potential(vec):
    return sp.expand(LAMBDA * (norm2(vec) - V**2/2)**2)


def edge_energy(phi,links,edge):
    i,j=edge
    return norm2(phi[i]-links[(i,j)]*phi[j])


def action(phi,links):
    assert len(phi) == N
    return sp.simplify(KAPPA/A * sum((edge_energy(phi,links,e) for e in EDGES),sp.Integer(0))
                       + A*sum((potential(p) for p in phi),sp.Integer(0)))


def local_gauge(phi,links,transform):
    """Combined SU(2)_L x U(1)_Y transformation in Higgs representation."""
    assert len(transform) == N and len(phi) == N
    out_phi = tuple(transform[i]*phi[i] for i in range(N))
    out_links = {(i,j):transform[i]*links[(i,j)]*transform[j].H for i,j in EDGES}
    return out_phi, out_links


def exact_site_gauges():
    """Predetermined 7 local group elements, not adjusted to the seven fields."""
    # SU(2) elements with rational entries and local U(1) phases.
    rot = [I2,R_SU2,R_SU2.T,sp.Matrix([[0,-1],[1,0]]),
           sp.Matrix([[0,1],[-1,0]]),-I2,sp.Matrix([[1,0],[0,-1]])*sp.I]
    # The last element diag(i,-i) is SU(2), as det=1.
    phase = [sp.Integer(1),Q_U1,sp.conjugate(Q_U1),-1,sp.I,-sp.I,sp.Integer(1)]
    H = tuple(phase[i]*rot[i] for i in range(N))
    assert all(sp.simplify((h.H*h-I2))==sp.zeros(2) for h in H)
    assert all(sp.simplify(rot[i].det()-1)==0 for i in range(N))
    return H


def wilson_trace(links):
    M=I2
    for e in EDGES:
        M=M*links[e]
    return sp.simplify(sp.trace(M))


def lattice_distance_energy_from_K(K):
    return sum((sp.Rational(K[i][j])/4 for i,j in EDGES),sp.Integer(0))


def case_distance_matrices():
    return {
        'Polybrutus_3D': geom.matrix_from_coordinates(geom.POLY3,geom.W3),
        '4D_A': geom.matrix_from_coordinates(geom.Q4_A,geom.W4_A),
        '4D_B_not_on_vacuum_sphere': geom.K4_B,
        '4D_C': geom.K4_C,
        'Cube_plus_1': geom.matrix_from_coordinates(geom.CUBE7,geom.W_CUBE),
    }


def orientations():
    return gauge.input_cases()


def local_energy_analytic(phi):
    """Exact result from K and identity links, avoiding large radicals."""
    assert len(phi)==N
    return sp.simplify(sum((norm2(phi[i]-phi[j]) for i,j in EDGES),sp.Integer(0)))


def prescribed_holonomy_energy(phi):
    """Predeclared link twist on the final edge, without fitting to any field."""
    i,j=EDGES[-1]
    return sp.simplify(local_energy_analytic(phi)-norm2(phi[i]-phi[j])
                       +norm2(phi[i]-TWIST*phi[j]))


def stress_control():
    """No mass claim: local V(h) and its curvature at the illustrative vacuum."""
    h=sp.symbols('h',real=True)
    q=(V+h)/sp.sqrt(2)
    pot=sp.expand(LAMBDA * (q*q - V**2/2)**2)
    expected=sp.expand(LAMBDA * (2*V*h+h*h)**2/4)
    assert sp.simplify(pot-expected)==0
    assert pot.subs(h,0)==0
    m2=sp.simplify(sp.diff(pot,h,2).subs(h,0))
    assert m2 == 2*LAMBDA*V**2 == 4
    return {'onsite_potential':str(pot),'local_curvature_m2_toy_units':str(m2)}


def as_text_number(value):
    return str(sp.simplify(value))


def run():
    check_fixed_protocol()
    # Reuse independent exact source certificates; 4D B is not on S3.
    prior,_ = geom.run()
    assert prior['4D_B']['spherical'] is False
    matrices=case_distance_matrices()
    spinors=orientations()
    names={'Polybrutus_3D':'Polybrutus_3D_embedded_S3',
           '4D_A':'4D_A_on_S3','4D_C':'4D_C_on_S3',
           'Cube_plus_1':'Cube_six_faces_plus_one_S3'}
    energies={}
    for name,spin_key in names.items():
        phi=spinors[spin_key]
        assert len(phi)==N and all(sp.simplify(norm2(p)-1)==0 for p in phi)
        if name!='Cube_plus_1':
            assert geom.histogram(matrices[name])==geom.EXPECTED
        fixed = lattice_distance_energy_from_K(matrices[name])
        direct = local_energy_analytic(phi)
        assert sp.simplify(fixed-direct)==0
        twist = prescribed_holonomy_energy(phi)
        assert twist.is_real and sp.simplify(twist)>=0
        energies[name]={
            'identity_link_action':as_text_number(fixed),
            'predeclared_twisted_link_action':as_text_number(twist),
            'potential_at_vacuum':'0',
            'cycle_edges':N,
        }
    # B is admitted as geometric witness but cannot be a 7-site vacuum state.
    energies['4D_B_not_on_vacuum_sphere']={
        'identity_link_distance_sum_only':str(lattice_distance_energy_from_K(matrices['4D_B_not_on_vacuum_sphere'])),
        'included_in_fixed_norm_Higgs_comparison':False,
    }
    phi=spinors['Cube_six_faces_plus_one_S3']
    H=exact_site_gauges()
    W=unit_links(twisted=True)
    new_phi,new_W=local_gauge(phi,W,H)
    for edge in EDGES:
        assert sp.simplify(edge_energy(phi,W,edge)-edge_energy(new_phi,new_W,edge))==0
    assert sp.simplify(action(phi,W)-action(new_phi,new_W))==0
    assert sp.simplify(wilson_trace(W)-wilson_trace(new_W))==0
    assert wilson_trace(unit_links())==2
    assert wilson_trace(W)==2*sp.Rational(3,5)*Q_U1
    # Explicit algebraic 7-link gauge covariance across EVERY cospherical case.
    for case_name, spin_key in names.items():
        state = spinors[spin_key]
        p2,w2 = local_gauge(state,W,H)
        for edge in EDGES:
            assert sp.simplify(edge_energy(state,W,edge)
                               - edge_energy(p2,w2,edge)) == 0, (case_name,edge)
        assert sp.simplify(action(state,W)-action(p2,w2)) == 0, case_name
        assert sp.simplify(wilson_trace(W)-wilson_trace(w2)) == 0, case_name
    # Control: the assignment of field states to sites is NOT specified by
    # the unlabeled 21 distances. Swapping 3 and 4 changes cycle energy.
    source=matrices['Polybrutus_3D']
    assignment=list(range(N))
    assignment[3],assignment[4] = assignment[4],assignment[3]
    Kswap=[[source[assignment[i]][assignment[j]] for j in range(N)]
           for i in range(N)]
    original_cycle=lattice_distance_energy_from_K(source)
    swapped_cycle=lattice_distance_energy_from_K(Kswap)
    assert original_cycle == sp.Rational(31,2)
    assert swapped_cycle == sp.Rational(37,2)
    assert Counter(source[i][j] for i,j in combinations(range(N),2)) == \
           Counter(Kswap[i][j] for i,j in combinations(range(N),2))
    # Incorrect control: transform the field and leave links unchanged.
    changed_only=sp.simplify(action(new_phi,W))
    assert sp.simplify(changed_only-action(phi,W)) != 0
    # Pure-gauge connections fitted to phi make all covariant gradients zero.
    # This demonstrates non-identifiability, NOT a competing physical model.
    frames=tuple(gauge.su2_from_spinor(p) for p in phi)
    pure={(i,j):frames[i]*frames[j].H for i,j in EDGES}
    assert sp.simplify(action(phi,pure))==0
    return {
        'status':'PASS',
        'model':'1D closed seven-site SU2xU1Y static toy; fixed site map',
        'parameters':{'N':7,'a':'1','kappa':'1','lambda':'1','v':'sqrt(2)',
                      'g_SU2':'1','g_U1':'1','Y_Higgs':'1/2',
                      'background_SU2':'[[3/5,-4/5],[4/5,3/5]]',
                      'background_U1':'(3+4*i)/5','twist_edge':'6 -> 0'},
        'cases':energies,
        'wilson_trace_identity':'2',
        'wilson_trace_twisted':str(wilson_trace(W)),
        'gauge_covariance_check':'4 of 4 cases: 7 edges each + full action + Wilson trace',
        'label_permutation_negative_control':{
            'original_cycle_energy':'31/2','after_swap_sites_3_4':'37/2',
            'same_21_unlabeled_distances':True},
        'wrong_field_only_change_detected':True,
        'pure_gauge_vacuum_energy':'0',
        'onsite':stress_control(),
        'interpretation':'toy geometry/gauge test; no particle or Higgs mass predicted',
    }


if __name__=='__main__':
    print(json.dumps(run(),indent=2,ensure_ascii=False))
