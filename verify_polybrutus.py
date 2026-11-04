#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Polybrutus: vérificateur indépendant exact, Python standard seulement.

Coordonnées normalisées: x = a*sqrt(3), y = b, z = c
avec a, b, c dans Q (fractions exactes). Les carrés de distances
sont alors des rationnels exacts, sans calcul flottant.
"""

from collections import Counter
from fractions import Fraction as F
from itertools import combinations
import argparse
import json

# Noms et triplets (coefficient de sqrt(3) sur X, Y, Z), R=1.
POINTS = (
    ("P1 Nord", F(0), F(0), F(1)),
    ("P2 Haut-1", F(1, 2), F(0), F(1, 2)),
    ("P3 Haut-2", F(-1, 4), F(3, 4), F(1, 2)),
    ("P4 Haut-3", F(-1, 4), F(-3, 4), F(1, 2)),
    ("P5 Sud", F(0), F(0), F(-1)),
    ("P6 Bas-1", F(1, 2), F(0), F(-1, 2)),
    ("P7 Bas-2", F(-1, 2), F(0), F(-1, 2)),
)
EXPECTED_MULTIPLICITIES = {4: 6, 7: 2, 9: 3, 12: 6, 13: 2, 16: 2}

def squared_norm(p):
    return 3 * p[1] ** 2 + p[2] ** 2 + p[3] ** 2

def squared_distance(p, q):
    return 3 * (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2 + (p[3] - q[3]) ** 2

def verify():
    assert len(POINTS) == 7, "Il faut exactement sept sommets"
    assert len({tuple(p[1:]) for p in POINTS}) == 7, "Sommets distincts exigés"
    assert all(squared_norm(p) == 1 for p in POINTS), "Rayons incorrects"

    edges = list(combinations(POINTS, 2))
    assert len(edges) == 21, "Nombre de liaisons incorrect"
    scaled_k = [4 * squared_distance(p, q) for p, q in edges]
    assert all(k.denominator == 1 for k in scaled_k), "k non entier"
    multiplicities = dict(sorted(Counter(int(k) for k in scaled_k).items()))
    assert multiplicities == EXPECTED_MULTIPLICITIES, "Familles de distances incorrectes"

    barycentric = (
        sum((p[1] for p in POINTS), F(0)) / 7,
        sum((p[2] for p in POINTS), F(0)) / 7,
        sum((p[3] for p in POINTS), F(0)) / 7,
    )
    assert barycentric == (F(0), F(0), F(1, 14)), "Barycentre incorrect"

    assert all(p[3] == F(1, 2) for p in POINTS[1:4]), "Niveau supérieur incorrect"
    assert all(p[3] == F(-1, 2) for p in POINTS[5:7]), "Niveau inférieur incorrect"
    assert POINTS[5][2] == POINTS[6][2] == 0, "Les deux points bas doivent être dans XZ"
    assert POINTS[0][3] == 1 and POINTS[4][3] == -1, "Pôles incorrects"

    # À z = ±1/2, r_horizontal²=3/4, donc cos(angle) = 1/2 : angle = 60°.
    assert all(3*p[1]**2 + p[2]**2 == F(3,4)
               for p in (*POINTS[1:4], *POINTS[5:7])), "Angle 60° incorrect"

    return {
        "status": "PASS",
        "sommets": len(POINTS),
        "liaisons": len(edges),
        "sur_sphere_unitaire": True,
        "angle_polaire": "60 degrés",
        "k_multiplicites": multiplicities,
        "barycentre_normalise": ["0", "0", "1/14"],
    }

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="Sortie JSON reproductible")
    args = parser.parse_args()
    result = verify()
    if args.json:
        print(json.dumps(result, ensure_ascii=False, indent=2))
    else:
        print("POLYBRUTUS : PASS (arithmétique exacte, sans dépendances)")
        print("7 sommets | 21 liaisons | rayon constant | angle polaire 60°")
        print("k -> multiplicité:", result["k_multiplicites"])
        print("Barycentre normalisé : (0, 0, 1/14)")

if __name__ == "__main__":
    main()
