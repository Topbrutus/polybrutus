#!/usr/bin/env python3
# coding: utf-8
import itertools
import math
import unittest
from collections import Counter
from fractions import Fraction
from pathlib import Path
import numpy as np
import sympy as sp
import electrofaible_sept_sites as lab

class TestElectrofaibleSeptSites(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cases=lab.load_spinors()

    def test_01_inputs_present(self):
        self.assertEqual(set(self.cases),{'Polybrutus_3D_embedded_S3','4D_A_on_S3','4D_C_on_S3','Cube_six_faces_plus_one_S3'})

    def test_02_unit_spinors(self):
        for p in self.cases.values():
            self.assertEqual(p.shape,(7,2))
            self.assertTrue(np.allclose(np.sum(abs(p)**2,axis=1),1))

    def test_03_exact_multisets(self):
        import json
        raw=json.loads(Path(lab.__file__).with_name('GEOMETRIES_EXACTES.json').read_text(encoding='utf8'))['spinors']
        expected=Counter({sp.Integer(1):6,sp.Rational(7,4):2,sp.Rational(9,4):3,
                          sp.Integer(3):6,sp.Rational(13,4):2,sp.Integer(4):2})
        for name in ('Polybrutus_3D_embedded_S3','4D_A_on_S3','4D_C_on_S3'):
            p=[sp.Matrix(list(map(sp.sympify,row))) for row in raw[name]]
            h=Counter(sp.simplify(((p[i]-p[j]).conjugate().T*(p[i]-p[j]))[0])
                      for i,j in itertools.combinations(range(7),2))
            self.assertEqual(h,expected,name)

    def test_04_exact_vacuum(self):
        import json
        raw=json.loads(Path(lab.__file__).with_name('GEOMETRIES_EXACTES.json').read_text(encoding='utf8'))['spinors']
        for name,rows in raw.items():
            for row in rows:
                self.assertEqual(sp.simplify(sum(sp.conjugate(sp.sympify(z))*sp.sympify(z) for z in row)),1,name)

    def test_05_background_su2_det_and_hyper_unit(self):
        for bg in lab.BACKGROUNDS:
            S,Z=lab.links(bg)
            for e in lab.EDGES:
                self.assertTrue(np.allclose(S[e].conj().T@S[e],lab.I2))
                self.assertAlmostEqual(abs(np.linalg.det(S[e])),1)
                self.assertAlmostEqual(abs(Z[e]),1)

    def test_06_loop_identity(self):
        S,Z=lab.links('identity');su,hy,combined=lab.wilson(S,Z)
        self.assertTrue(np.allclose(su,lab.I2))
        self.assertEqual(hy,1)
        self.assertTrue(np.allclose(combined,lab.I2))

    def test_07_su2_wilson_nontrivial(self):
        S,Z=lab.links('su2_only');su,hy,combined=lab.wilson(S,Z)
        self.assertAlmostEqual(np.trace(su).real,6/5)
        self.assertEqual(hy,1)
        self.assertTrue(np.allclose(combined,lab.SQ))

    def test_08_neutral_em_both_holonomies_nontrivial(self):
        S,Z=lab.links('neutral_em');su,hy,combined=lab.wilson(S,Z)
        self.assertFalse(np.allclose(su,lab.I2))
        self.assertNotEqual(hy,1)
        self.assertTrue(np.allclose(combined@lab.E2,lab.E2))
        self.assertTrue(np.allclose(combined,np.diag([lab.Q**2,1])))

    def test_09_neutral_vev_zero_action(self):
        S,Z=lab.links('neutral_em');p=np.tile(lab.E2,(7,1))
        self.assertAlmostEqual(lab.energy(p,S,Z),0)

    def test_10_charged_state_fails_zero(self):
        S,Z=lab.links('neutral_em');p=np.tile(lab.E1,(7,1))
        self.assertAlmostEqual(lab.energy(p,S,Z),64/25)

    def test_11_generic_vev_positive(self):
        S,Z=lab.links('su2_only');p=np.tile(lab.E2,(7,1))
        self.assertAlmostEqual(lab.energy(p,S,Z),4/5)

    def test_12_generic_rigorous_lower_bound(self):
        S,Z=lab.links('su2_only')
        W=lab.wilson(S,Z)[2]
        self.assertTrue(np.allclose((lab.I2-W).conj().T@(lab.I2-W),4/5*lab.I2))
        for p in self.cases.values():
            self.assertGreaterEqual(lab.energy(p,S,Z),4/35-1e-11)

    def test_13_operator_hermitian(self):
        for bg in lab.BACKGROUNDS:
            S,Z=lab.links(bg);L=lab.twisted_laplacian(S,Z)
            self.assertEqual(L.shape,(14,14))
            self.assertTrue(np.allclose(L,L.conj().T))
            self.assertGreaterEqual(np.linalg.eigvalsh(L)[0],-1e-12)

    def test_14_operator_energy_identity(self):
        for bg in lab.BACKGROUNDS:
            S,Z=lab.links(bg);L=lab.twisted_laplacian(S,Z)
            for p in self.cases.values():
                self.assertAlmostEqual(lab.energy(p,S,Z),float(np.vdot(p.reshape(14),L@p.reshape(14)).real),places=10)

    def test_15_spectrum_exact_minimum(self):
        for bg in lab.BACKGROUNDS:
            S,Z=lab.links(bg);eig=np.linalg.eigvalsh(lab.twisted_laplacian(S,Z))
            self.assertAlmostEqual(7*eig[0],lab.minimum_energy_analytical(bg),places=10)

    def test_16_explicit_minimizer(self):
        p=lab.generic_minimizer();S,Z=lab.links('su2_only')
        self.assertTrue(np.allclose(np.sum(abs(p)**2,axis=1),1))
        self.assertAlmostEqual(lab.energy(p,S,Z),lab.minimum_energy_analytical('su2_only'),places=11)

    def test_17_local_full_electroweak_gauge(self):
        su,hy=lab.local_gauges()
        for p in self.cases.values():
            for bg in lab.BACKGROUNDS:
                S,Z=lab.links(bg)
                p2,S2,Z2=lab.local_transform(p,S,Z,su,hy)
                self.assertAlmostEqual(lab.energy(p,S,Z),lab.energy(p2,S2,Z2),places=10)
                old=lab.wilson(S,Z);new=lab.wilson(S2,Z2)
                self.assertTrue(np.allclose(new[0],su[0]@old[0]@su[0].conj().T))
                self.assertTrue(np.allclose(new[1],old[1]))
                self.assertAlmostEqual(np.trace(new[2]).real,np.trace(old[2]).real)
                self.assertAlmostEqual(np.trace(new[2]).imag,np.trace(old[2]).imag)

    def test_18_wrong_field_only_control(self):
        p=self.cases['Polybrutus_3D_embedded_S3'];S,Z=lab.links('neutral_em');su,hy=lab.local_gauges()
        changed,_,_=lab.local_transform(p,S,Z,su,hy)
        self.assertGreater(abs(lab.energy(changed,S,Z)-lab.energy(p,S,Z)),1e-4)

    def test_19_relax_connections_erases_poly_geometry(self):
        for p in self.cases.values():
            S,Z=lab.fitted_flat_links(p)
            self.assertLess(lab.energy(p,S,Z),1e-23)
            self.assertTrue(np.allclose(lab.wilson(S,Z)[2],lab.I2))

    def test_20_relaxed_pure_gauge_is_original_vacuum(self):
        for p in self.cases.values():
            S,Z=lab.fitted_flat_links(p)
            frames=lab.from_neutral(p)
            su=np.array([f.conj().T for f in frames]);hy=np.ones(7,dtype=complex)
            pp,SS,ZZ=lab.local_transform(p,S,Z,su,hy)
            self.assertTrue(np.allclose(pp,np.tile(lab.E2,(7,1))))
            self.assertTrue(all(np.allclose(SS[e],lab.I2) and abs(ZZ[e]-1)<1e-12 for e in lab.EDGES))

    def test_21_prior_cycle_benchmark(self):
        S,Z=lab.links('identity')
        expected={'Polybrutus_3D_embedded_S3':31/2,'4D_A_on_S3':18,'4D_C_on_S3':16,'Cube_six_faces_plus_one_S3':14}
        for name,v in expected.items():
            self.assertAlmostEqual(lab.energy(self.cases[name],S,Z),v,places=10)

    def test_22_all_permutation_means(self):
        for name in ('Polybrutus_3D_embedded_S3','4D_A_on_S3','4D_C_on_S3'):
            d=lab.identity_cycle_energies(self.cases[name])
            self.assertEqual(d['count'],5040)
            self.assertAlmostEqual(d['mean'],65/4,places=9)
            self.assertAlmostEqual(d['predicted_mean'],65/4,places=9)
        self.assertNotEqual(lab.identity_cycle_energies(self.cases['Cube_six_faces_plus_one_S3'])['predicted_mean'],65/4)

    def test_23_relabel_changes_energy(self):
        p=self.cases['Polybrutus_3D_embedded_S3'];S,Z=lab.links('identity')
        pp=p.copy();pp[[3,4]]=pp[[4,3]]
        self.assertAlmostEqual(lab.energy(p,S,Z),31/2)
        self.assertAlmostEqual(lab.energy(pp,S,Z),37/2)

    def test_24_exact_broken_su2_analytical_bound(self):
        q=sp.Rational(3,5)+sp.I*sp.Rational(4,5)
        R=sp.diag(q,sp.conjugate(q));W=q*R
        self.assertEqual(sp.simplify(R.det()),1)
        self.assertEqual(sp.simplify((R.H*R-sp.eye(2))),sp.zeros(2))
        self.assertEqual(sp.simplify(W*sp.Matrix([0,1])),sp.Matrix([0,1]))
        self.assertEqual(sp.simplify((sp.eye(2)-R).H*(sp.eye(2)-R)),sp.Rational(4,5)*sp.eye(2))
        self.assertEqual(sp.simplify(((sp.eye(2)-W)*sp.Matrix([1,0])).norm()**2),sp.Rational(64,25))

    def test_25_predeclared_inputs_and_json(self):
        r=lab.run()
        self.assertEqual(r['status'],'PASS')
        self.assertEqual(len(r['results']),4)
        self.assertEqual(len(r['backgrounds']),4)

if __name__=='__main__':
    unittest.main()
