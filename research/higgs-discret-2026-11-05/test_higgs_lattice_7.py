#!/usr/bin/env python3
# Independent unit test suite for the fixed seven-site gauge model.
import unittest
from collections import Counter
from itertools import combinations
import sympy as s
import higgs_lattice_7 as m


class DiscreteHiggsSevenSiteTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.results=m.run()  # All four complete certificates and negative controls.
        cls.state=m.orientations()
        cls.K=m.case_distance_matrices()
        cls.H=m.exact_site_gauges()
        cls.W=m.unit_links(twisted=True)

    def test_01_predeclared_parameters(self):
        self.assertTrue(m.check_fixed_protocol())
        self.assertEqual(m.N,7)
        self.assertEqual(m.Y_HIGGS,s.Rational(1,2))

    def test_02_positions_are_spatial_and_separate(self):
        positions=m.positions()
        self.assertEqual(len(positions),7)
        self.assertEqual(len(set(positions)),7)
        self.assertEqual(s.simplify(2*m.RADIUS*s.sin(s.pi/7)),1)

    def test_03_graph_is_one_closed_cycle(self):
        self.assertEqual(m.EDGES,((0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(6,0)))

    def test_04_cospherical_input_norms(self):
        for k,phi in self.state.items():
            self.assertEqual(len(phi),7,k)
            for p in phi:
                self.assertEqual(s.simplify(m.norm2(p)),1,k)

    def test_05_exact_histogram_preserved(self):
        for name in ('Polybrutus_3D','4D_A','4D_B_not_on_vacuum_sphere','4D_C'):
            self.assertEqual(m.geom.histogram(self.K[name]),m.geom.EXPECTED)

    def test_06_21_pairwise_sums_identical(self):
        sums={sum((s.Rational(K[i][j],4) for i,j in combinations(range(7),2)),s.Integer(0))
              for name,K in self.K.items() if name!='Cube_plus_1'}
        self.assertEqual(sums,{s.Rational(195,4)})

    def test_07_cycle_action_differs_by_class(self):
        cases=self.results['cases']
        self.assertEqual([cases[n]['identity_link_action'] for n in
                          ('Polybrutus_3D','4D_A','4D_C','Cube_plus_1')],
                         ['31/2','18','16','14'])

    def test_08_fixed_twist_not_equal_identity(self):
        cases=self.results['cases']
        for n in ('Polybrutus_3D','4D_A','4D_C','Cube_plus_1'):
            self.assertNotEqual(cases[n]['identity_link_action'],
                                cases[n]['predeclared_twisted_link_action'])

    def test_09_site_permutation_negative_control(self):
        ctr=self.results['label_permutation_negative_control']
        self.assertEqual(ctr['original_cycle_energy'],'31/2')
        self.assertEqual(ctr['after_swap_sites_3_4'],'37/2')
        self.assertTrue(ctr['same_21_unlabeled_distances'])

    def test_10_su2_and_u1_composite_gauge_unitary(self):
        for h in self.H:
            self.assertEqual((h.H*h-m.I2).applyfunc(s.simplify),s.zeros(2))
        self.assertEqual(s.det(m.R_SU2),1)
        self.assertEqual(s.expand(m.Q_U1*s.conjugate(m.Q_U1)),1)

    def test_11_wilson_loop_is_gauge_invariant(self):
        initial=m.wilson_trace(self.W)
        self.assertEqual(initial,s.Rational(18,25)+s.Rational(24,25)*s.I)
        for phi in self.state.values():
            _,altered=m.local_gauge(phi,self.W,self.H)
            self.assertEqual(s.simplify(m.wilson_trace(altered)-initial),0)

    def test_12_seven_covariant_links_for_all_four_cases(self):
        for name,phi in self.state.items():
            new_phi,new_W=m.local_gauge(phi,self.W,self.H)
            for edge in m.EDGES:
                self.assertEqual(s.simplify(m.edge_energy(phi,self.W,edge)
                                            -m.edge_energy(new_phi,new_W,edge)),0,(name,edge))

    def test_13_complete_action_invariance(self):
        for name,phi in self.state.items():
            new_phi,new_W=m.local_gauge(phi,self.W,self.H)
            self.assertEqual(s.simplify(m.action(phi,self.W)-m.action(new_phi,new_W)),0,name)

    def test_14_onsite_potential_vanishes_for_all_vacua(self):
        for name,phi in self.state.items():
            for p in phi:
                self.assertEqual(s.simplify(m.potential(p)),0,name)

    def test_15_negative_field_only_change_breaks_invariance(self):
        phi=self.state['Cube_six_faces_plus_one_S3']
        new_phi,_=m.local_gauge(phi,self.W,self.H)
        self.assertNotEqual(s.simplify(m.action(phi,self.W)-m.action(new_phi,self.W)),0)

    def test_16_pure_gauge_fits_any_normalized_orientation(self):
        phi=self.state['Polybrutus_3D_embedded_S3']
        frames=[m.gauge.su2_from_spinor(p) for p in phi]
        links={(i,j):frames[i]*frames[j].H for i,j in m.EDGES}
        self.assertEqual(s.simplify(m.action(phi,links)),0)

    def test_17_local_potential_curvature_no_mass_prediction(self):
        self.assertEqual(m.stress_control()['local_curvature_m2_toy_units'],'4')
        self.assertIn('no particle',self.results['interpretation'])

    def test_18_noncospherical_control_not_claimed_higgs_vacuum(self):
        self.assertFalse(self.results['cases']['4D_B_not_on_vacuum_sphere']
                         ['included_in_fixed_norm_Higgs_comparison'])

    def test_19_wilson_distinguishes_fixed_backgrounds(self):
        self.assertEqual(m.wilson_trace(m.unit_links()),2)
        self.assertNotEqual(m.wilson_trace(self.W),2)


if __name__=='__main__':
    unittest.main(verbosity=2)
