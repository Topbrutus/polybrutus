#!/usr/bin/env python3
"""Minimal exact coordinate reference for Polybrutus v0.1.0 (unit radius).

This module exposes only the geometry required by the 3D/4D reconstruction
certificates. It deliberately omits unrelated experimental spectral/X7 data.
"""
from fractions import Fraction as F

# (x/sqrt(3), y, z) in exact rational coordinates.
EXACT_POINTS = (
    (F(0), F(0), F(1)),
    (F(1, 2), F(0), F(1, 2)),
    (F(-1, 4), F(3, 4), F(1, 2)),
    (F(-1, 4), F(-3, 4), F(1, 2)),
    (F(0), F(0), F(-1)),
    (F(1, 2), F(0), F(-1, 2)),
    (F(-1, 2), F(0), F(-1, 2)),
)

def d2e(i, j):
    """Exact squared distance between reference points i and j."""
    p, q = EXACT_POINTS[i], EXACT_POINTS[j]
    return 3 * (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2
