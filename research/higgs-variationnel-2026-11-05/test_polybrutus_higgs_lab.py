#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Reproducible unit tests for exact Polybrutus–Higgs exploration."""
import unittest
from collections import Counter
from fractions import Fraction as F
from itertools import combinations

import polybrutus_higgs_lab as lab


class TestPolybrutusHiggsLab(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cases, cls.control = lab.run()
        cls.D = {
            'Polybrutus_3D': lab.matrix_from_coordinates(lab.POLY3, lab.W3),
            '4D_A': lab.matrix_from_coordinates(lab.Q4_A, lab.W4_A),
            '4D_B': lab.K4_B,
            '4D_C': lab.K4_C,
        }

    def test_same_21_unlabeled_distances(self):
        for name, K in self.D.items():
            with self.subTest(name=name):
                self.assertEqual(lab.histogram(K), lab.EXPECTED)
                self.assertEqual(sum(K[i][j] for i,j in combinations(range(7),2)), 195)

    def test_affine_ranks(self):
        self.assertEqual([d['rank'] for d in self.cases.values()], [3,4,4,4])
        self.assertEqual([lab.rank(lab.centered_gram(K)) for K in self.D.values()], [3,4,4,4])

    def test_sphere_restriction(self):
        self.assertEqual([d['spherical'] for d in self.cases.values()], [True,True,False,True])
        for name in ('Polybrutus_3D','4D_A','4D_C'):
            self.assertEqual(self.cases[name]['sphere_radius_sq_anchor'], 1)

    def test_pairwise_sum_cannot_select(self):
        for name,d in self.cases.items():
            with self.subTest(name=name):
                self.assertEqual(d['E2'],F(6329,546))
                self.assertEqual(lab.pair_energy(self.D[name],2),
                                 sum((F(4,k)**2 * n for k,n in lab.EXPECTED.items()), F(0)))

    def test_triple_descriptor_distinguishes_all_four(self):
        assert_vals = [F(24435,64), F(23835,64), F(26817,64), F(24273,64)]
        actual = [d['E3'] for d in self.cases.values()]
        self.assertEqual(actual, assert_vals)
        self.assertEqual(len(set(actual)),4)

    def test_triangle_structures_are_not_equal(self):
        self.assertEqual(self.cases['Polybrutus_3D']['equi'],
                         Counter({F(3):1,F(9,4):1}))
        self.assertEqual(self.cases['4D_A']['equi'], Counter({F(1):1}))
        self.assertEqual(self.cases['4D_B']['equi'], Counter({F(1):1}))
        self.assertEqual(self.cases['4D_C']['equi'], Counter({F(3):1}))

    def test_no_stationary_simple_repulsion(self):
        for name in ('Polybrutus_3D','4D_A','4D_C'):
            self.assertGreater(self.cases[name]['projected_force_norm_sq'],0)
        self.assertIsNone(self.cases['4D_B']['projected_force_norm_sq'])

    def test_cube_control_not_polybrutus(self):
        cube = lab.matrix_from_coordinates(lab.CUBE7,lab.W_CUBE)
        self.assertEqual(lab.histogram(cube),Counter({F(4):6,F(6):12,F(12):3}))
        self.assertEqual(sum(F(c[0],7) for c in lab.CUBE7),F(4,7))
        self.assertEqual(self.control['cube_E2'],15)
        self.assertEqual(self.control['cube_E3'],135)
        self.assertNotEqual(lab.histogram(cube),lab.EXPECTED)

    def test_radial_higgs_potential_algebra(self):
        # λ=v=1 -> V(h)=1/4*(2h+h²)^2
        # Polynomial coefficients h^0,h^1,h^2,h^3,h^4 are (0,0,1,1,1/4).
        coeffs={0:F(0),1:F(0),2:F(1),3:F(1),4:F(1,4)}
        self.assertEqual(coeffs[0],0)
        self.assertEqual(coeffs[1],0)
        self.assertEqual(2*coeffs[2],2) # m_h^2/(λv²)
        self.assertEqual(self.control['higgs_radial_mass_squared_over_lambda_v_squared'],2)
        self.assertEqual(self.control['higgs_angular_modes_with_zero_potential_curvature'],3)

    def test_positive_third_moment_discriminant(self):
        self.assertEqual([d['moment3'] for d in self.cases.values()], [
            F(4629285,87808),F(5234925,87808),
            F(3972321,87808),F(4721139,87808)])


if __name__=='__main__':
    unittest.main(verbosity=2)
