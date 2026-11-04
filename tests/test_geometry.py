"""Tests unitaires des invariants du Polybrutus."""
import sys
from pathlib import Path
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from verify_polybrutus import POINTS, squared_norm, squared_distance, verify

class TestPolybrutus(unittest.TestCase):
    def test_unicity(self):
        self.assertEqual(len(POINTS), 7)
        self.assertEqual(len(set(tuple(p[1:]) for p in POINTS)), 7)

    def test_same_sphere(self):
        self.assertTrue(all(squared_norm(p) == 1 for p in POINTS))

    def test_six_exact_families(self):
        self.assertEqual(verify()["k_multiplicites"],
            {4:6, 7:2, 9:3, 12:6, 13:2, 16:2})

    def test_barycenter(self):
        self.assertEqual(verify()["barycentre_normalise"], ["0","0","1/14"])

    def test_polar_angle(self):
        self.assertEqual(verify()["angle_polaire"], "60 degrés")

    def test_poles(self):
        self.assertEqual(squared_distance(POINTS[0], POINTS[4]), 4)

if __name__ == "__main__":
    unittest.main()
