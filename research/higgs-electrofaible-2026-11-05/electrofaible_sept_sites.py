#!/usr/bin/env python3
# coding: utf-8
"""Sept-site *static* SU(2)xU(1) Higgs toy with a fixed external gauge background.

The experiment is deliberately one spatial dimension, dimensionless and NOT
Lorentz-invariant (nor a physical Standard Model simulation). A nonzero
Higgs mass or a new particle is not inferred from these computations.
"""
from __future__ import annotations
import itertools
import json
import math
from collections import Counter
from pathlib import Path
import numpy as np
import sympy as sp

N = 7
EDGES = tuple((i, (i+1) % N) for i in range(N))
Q = (3 + 4j) / 5  # independently fixed, exact rational components
SQ = np.diag([Q, Q.conjugate()]).astype(complex)
I2 = np.eye(2, dtype=complex)
E1 = np.array([1,0], dtype=complex)  # electrically charged upper state
E2 = np.array([0,1], dtype=complex)  # neutral lower Higgs direction, Q=T3+Y=0
LAM = 1.0  # not a physical lambda measurement
A = 1.0    # not a physical length in meters
KAPPA = 1.0

# Prespecified backgrounds. Each SU(2) holonomy is *the same* in
# su2_only and neutral_em; the independently prescribed U(1)Y changes.
BACKGROUNDS = ('identity', 'su2_only', 'neutral_em', 'hyper_only')

def load_spinors(path: Path | None = None) -> dict[str,np.ndarray]:
    path = path or Path(__file__).with_name('GEOMETRIES_EXACTES.json')
    src = json.loads(path.read_text(encoding='utf-8'))
    out = {}
    for name, rows in src['spinors'].items():
        phi = np.array([[complex(sp.N(sp.sympify(z), 18)) for z in row]
                        for row in rows], dtype=complex)
        if phi.shape != (N,2) or not np.allclose(np.linalg.norm(phi,axis=1),1,atol=1e-12):
            raise ValueError(f'Invalid S3 embedding for {name}')
        out[name] = phi
    return out


def links(background: str) -> tuple[dict,dict]:
    if background not in BACKGROUNDS:
        raise ValueError(background)
    S = {e: I2.copy() for e in EDGES}
    Z = {e: complex(1) for e in EDGES}
    if background in ('su2_only','neutral_em'):
        S[EDGES[-1]] = SQ.copy()
    if background in ('neutral_em','hyper_only'):
        Z[EDGES[-1]] = Q
    return S,Z


def energy(phi: np.ndarray, S: dict, Z: dict) -> float:
    if phi.shape != (N,2):
        raise ValueError('7 doublets C2 required')
    hop=sum(np.vdot(phi[i]-Z[e]*S[e]@phi[j],
                     phi[i]-Z[e]*S[e]@phi[j]).real for e in EDGES for i,j in (e,))
    pot=LAM*sum((np.vdot(v,v).real-1.0)**2 for v in phi)
    return float(KAPPA/A*hop + A*pot)


def wilson(S:dict,Z:dict) -> tuple[np.ndarray,complex,np.ndarray]:
    su=I2.copy(); hy=1+0j
    for e in EDGES:
        su=su@S[e]
        hy*=Z[e]
    return su,hy,hy*su


def local_gauges() -> tuple[np.ndarray,np.ndarray]:
    rot=np.array([[3/5,-4/5],[4/5,3/5]],complex)
    swap=np.array([[0,-1],[1,0]],complex)
    diag=SQ
    su=np.array([I2,rot,rot.conj().T,swap,swap.conj().T,diag,-I2])
    hy=np.array([1,Q,Q.conjugate(),-1,1j,-1j,Q**2],complex)
    assert all(np.allclose(h.conj().T@h,I2) and abs(np.linalg.det(h)-1)<1e-12 for h in su)
    assert np.allclose(np.abs(hy),1)
    return su,hy


def local_transform(phi:np.ndarray, S:dict, Z:dict,
                    su:np.ndarray,hy:np.ndarray) -> tuple[np.ndarray,dict,dict]:
    new_phi=np.stack([hy[i]*(su[i]@phi[i]) for i in range(N)])
    new_S={e:su[i]@S[e]@su[j].conj().T for e in EDGES for i,j in (e,)}
    new_Z={e:hy[i]*Z[e]*hy[j].conjugate() for e in EDGES for i,j in (e,)}
    return new_phi,new_S,new_Z


def from_neutral(phi:np.ndarray) -> np.ndarray:
    """SU2 frame F_i with F_i e2=phi_i, for each unit spinor."""
    frames=[]
    for a,b in phi:
        F=np.array([[b.conjugate(),a],[-a.conjugate(),b]],complex)
        assert np.allclose(F@E2,[a,b])
        assert np.allclose(F.conj().T@F,I2)
        assert abs(np.linalg.det(F)-1)<1e-11
        frames.append(F)
    return np.stack(frames)


def fitted_flat_links(phi:np.ndarray) -> tuple[dict,dict]:
    """CONTROL ONLY: choosing links *after seeing phi* erases any signal."""
    frames=from_neutral(phi)
    S={e: frames[i]@frames[j].conj().T for e in EDGES for i,j in (e,)}
    Z={e: 1+0j for e in EDGES}
    return S,Z


def twisted_laplacian(S:dict,Z:dict) -> np.ndarray:
    """14x14 Hermitian operator for the quadratic hopping functional."""
    L=np.zeros((2*N,2*N),complex)
    for (i,j) in EDGES:
        u=Z[(i,j)]*S[(i,j)]
        L[2*i:2*i+2,2*i:2*i+2]+=I2
        L[2*j:2*j+2,2*j:2*j+2]+=I2
        L[2*i:2*i+2,2*j:2*j+2]-=u
        L[2*j:2*j+2,2*i:2*i+2]-=u.conj().T
    assert np.allclose(L,L.conj().T)
    return L


def minimum_energy_analytical(background:str) -> float:
    """For seven unit doublets: the cyclic Laplacian min is attained by
    a uniform-amplitude eigenmode; the spectral lower bound is exact.
    """
    if background in ('identity','neutral_em'):return 0.0
    theta=math.atan2(4,3)
    return 2*N*(1-math.cos(theta/N))


def generic_minimizer() -> np.ndarray:
    """Uniform equal-phase-step minimizer for the SU2-only holonomy."""
    theta=math.atan2(4,3)
    return np.stack([np.exp(-1j*theta*i/N)*E2 for i in range(N)])


def identity_cycle_energies(phi:np.ndarray) -> dict:
    """Exhaustively sample all 5040 label-site assignments (already a prior control)."""
    dist=np.array([[np.vdot(phi[i]-phi[j],phi[i]-phi[j]).real for j in range(N)] for i in range(N)])
    values=[sum(dist[p[i],p[(i+1)%N]] for i in range(N)) for p in itertools.permutations(range(N))]
    distinct=len(set(round(v,10) for v in values))
    totalpairs=float(np.sum(np.triu(dist,1)))
    return {'count':len(values),'distinct_values':distinct,'minimum':min(values),
            'maximum':max(values),'mean':sum(values)/len(values),
            'predicted_mean':totalpairs/3}


def run() -> dict:
    cases=load_spinors(); results={}
    H,z=local_gauges()
    for name,phi in cases.items():
        calc={}
        for bg in BACKGROUNDS:
            S,Z=links(bg)
            before=energy(phi,S,Z)
            after_phi,after_S,after_Z=local_transform(phi,S,Z,H,z)
            after=energy(after_phi,after_S,after_Z)
            assert abs(before-after)<1e-10,(name,bg,before,after)
            w=wilson(S,Z);w2=wilson(after_S,after_Z)
            assert np.allclose(w2[2],z[0]*H[0]@w[2]@H[0].conj().T*z[0].conjugate())
            assert abs(np.trace(w[2])-np.trace(w2[2]))<1e-10
            calc[bg]=round(before,12)
        freeS,freeZ=fitted_flat_links(phi)
        freeE=energy(phi,freeS,freeZ)
        assert freeE<1e-22 and np.allclose(wilson(freeS,freeZ)[2],I2,atol=1e-11)
        results[name]={'fixed_assignment_energies':calc,
                       'relaxed_connection_energy':round(freeE,14),
                       'identity_cycle_permutations':identity_cycle_energies(phi)}
    spectral={}
    for bg in BACKGROUNDS:
        S,Z=links(bg)
        ev=np.linalg.eigvalsh(twisted_laplacian(S,Z))
        predicted=minimum_energy_analytical(bg)
        measured=N*float(ev[0])
        assert abs(predicted-measured)<1e-10,(bg,predicted,measured)
        spectral[bg]={
            'minimum_for_seven_unit_doublets':round(predicted,12),
            'lowest_laplacian_eigenvalue':round(float(ev[0]),12),
            'combined_wilson_trace':str(complex(round(np.trace(wilson(S,Z)[2]).real,12),round(np.trace(wilson(S,Z)[2]).imag,12))),
        }
    S,Z=links('neutral_em')
    assert energy(np.stack([E2]*N),S,Z)<1e-26
    assert abs(energy(np.stack([E1]*N),S,Z)-64/25)<1e-12
    S,Z=links('su2_only')
    assert abs(energy(np.stack([E2]*N),S,Z)-4/5)<1e-12
    assert abs(energy(generic_minimizer(),S,Z)-minimum_energy_analytical('su2_only'))<1e-12
    return {'status':'PASS','sites':N,'edges':N,'backgrounds':list(BACKGROUNDS),
            'independent_parameter_q':'(3+4i)/5',
            'model':'static 1D ring, SU(2)_L and U(1)_Y links in the Higgs representation, a=kappa=lambda=1, norm Phi=1',
            'results':results,'spectral':spectral,
            'neutral_vacuum_action':'0','charged_constant_action_neutral_em':'64/25',
            'generic_constant_action':'4/5',
            'generic_rigorous_lower_bound':'4/35',
            'generic_exact_minimum_formula':'14*(1-cos(atan2(4,3)/7))',
            'no_go':'all four unit-doublet input cases admit flat SU(2) links with E=0 when links are fitted to the fields',
            'interpretation':'no electroweak particle or Higgs mass claim; static gauge-invariant toy only'}

if __name__=='__main__':
    print(json.dumps(run(),ensure_ascii=False,indent=2))
