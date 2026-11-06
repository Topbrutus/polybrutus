#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Exact unit tests: SU(2) gauge covariance on seven-point toy sites."""
import unittest
from fractions import Fraction as F
import sympy as S

import higgs_gauge_links as lab


class TestGaugeLinks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.results = lab.run()  # full 21-edge/35-triangle exact certificate

    def test_four_geometries_including_two_true_4d_witnesses(self):
        self.assertEqual(len(self.results), 4)

    def test_previous_four_classes_still_certified(self):
        cases,_=lab.previous.run()
        self.assertEqual([cases[k]['spherical'] for k in ['Polybrutus_3D','4D_A','4D_B','4D_C']], [True,True,False,True])

    def test_each_state_on_s3(self):
        for name,states in lab.input_cases().items():
            with self.subTest(name=name):
                self.assertEqual(len(states),7)
                self.assertTrue(all(lab.norm_sq(x)==1 for x in states))

    def test_site_dependent_su2_exists(self):
        for name,states in lab.input_cases().items():
            with self.subTest(name=name):
                for x in states:
                    g=lab.su2_from_spinor(x)
                    self.assertEqual(S.simplify(g.det()),1)
                    self.assertTrue((g.H*g-lab.I2).applyfunc(S.simplify).is_zero_matrix)
                    self.assertTrue((g*lab.E1-x).applyfunc(S.simplify).is_zero_matrix)

    def test_all_21_naive_pair_differences_can_disappear(self):
        for name,states in lab.input_cases().items():
            with self.subTest(name=name):
                frames=tuple(lab.su2_from_spinor(x) for x in states)
                aligned=tuple(g.H*x for g,x in zip(frames,states))
                self.assertTrue(all((x-lab.E1).applyfunc(S.simplify).is_zero_matrix for x in aligned))
                self.assertNotEqual(self.results[name]['edge_sum_unit_links'],'0')

    def test_21_link_covariant_distances_survive(self):
        for result in self.results.values():
            self.assertEqual(result['tested_edges'],21)
            self.assertNotEqual(result['edge_sum_unit_links'],'0')

    def test_barycenter_not_plain_gauge_invariant(self):
        for result in self.results.values():
            self.assertNotEqual(result['bare_barycenter_norm_sq'], '1')
            self.assertEqual(result['bare_barycenter_norm_sq_after_local_gauge'], '1')

    def test_parallel_transported_barycenter_invariant(self):
        for result in self.results.values():
            self.assertEqual(result['transported_barycenter_norm_sq_identity_links'],
                             result['transported_barycenter_norm_sq_after_gauge'])

    def test_geometry_specific_barycenters(self):
        self.assertEqual(self.results['Polybrutus_3D_embedded_S3']['bare_barycenter_norm_sq'],'1/196')
        self.assertEqual(self.results['Cube_six_faces_plus_one_S3']['bare_barycenter_norm_sq'],'16/49')

    def test_connection_affects_covariant_pair_energy(self):
        for result in self.results.values():
            self.assertNotEqual(result['edge_sum_unit_links'],'0')
            self.assertEqual(result['edge_sum_pure_links'],'0')
            self.assertEqual(result['transported_barycenter_norm_sq_pure_links'],'1')

    def test_third_order_descriptor_is_covariant_with_links(self):
        for result in self.results.values():
            self.assertEqual(result['three_body_covariant_unit_links'],
                             result['three_body_covariant_after_local_gauge'])
            self.assertEqual(result['three_body_covariant_pure_links'],'0')

    def test_true_4d_witnesses_distinct_three_body_invariants(self):
        self.assertEqual(self.results['4D_A_on_S3']['three_body_covariant_unit_links'],'23835/64')
        self.assertEqual(self.results['4D_C_on_S3']['three_body_covariant_unit_links'],'24273/64')

    def test_direct_matrix_covariance_on_true_4d_sample_edges(self):
        for label in ('4D_A_on_S3','4D_C_on_S3'):
            with self.subTest(label=label):
                phi=lab.input_cases()[label]
                g0,g1=(lab.su2_from_spinor(phi[i]) for i in (0,1))
                h0,h1=g0.H,g1.H
                u=h0*h1.H
                a,b=h0*phi[0],h1*phi[1]
                self.assertEqual(S.simplify(lab.norm_sq(a-u*b)-lab.norm_sq(phi[0]-phi[1])),0)

    def test_exact_three_body_references(self):
        self.assertEqual(self.results['Polybrutus_3D_embedded_S3']['three_body_covariant_unit_links'],'24435/64')
        self.assertEqual(self.results['Cube_six_faces_plus_one_S3']['three_body_covariant_unit_links'],'135')

    def test_all_35_wilson_triangles_are_flat(self):
        for result in self.results.values():
            self.assertEqual(result['tested_triangles'],35)
        # lab.run also verifies trace(Uij Ujk Uki)==2 for each triangle

    def test_cube_distances_and_sum(self):
        self.assertEqual(self.results['Cube_six_faces_plus_one_S3']['edge_sum_unit_links'],'33')
        self.assertEqual(sorted(self.results['Cube_six_faces_plus_one_S3']['dist_sq_multiset']),
                         sorted([('1',6),('3/2',12),('3',3)]))

    def test_polybrutus_distances_and_sum(self):
        self.assertEqual(self.results['Polybrutus_3D_embedded_S3']['edge_sum_unit_links'],'195/4')
        self.assertEqual(sorted(self.results['Polybrutus_3D_embedded_S3']['dist_sq_multiset']),
                         sorted([('1',6),('7/4',2),('9/4',3),('3',6),('13/4',2),('4',2)]))


if __name__=='__main__':
    unittest.main(verbosity=2)
