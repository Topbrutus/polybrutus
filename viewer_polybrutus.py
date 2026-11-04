#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""POLYBRUTUS — Visualiseur 3D interactif (version indépendante).

Utilisation:
    pip install matplotlib numpy
    python viewer_polybrutus.py --radius 10
    python viewer_polybrutus.py --radius 10 --save polybrutus.png

Glisser avec la souris pour tourner la vue ; molette pour zoomer.
Aucun accès réseau ; aucune collecte de données.
"""
from argparse import ArgumentParser
from itertools import combinations
from math import pi, sqrt
from pathlib import Path

def coordinates(R):
    """Points analytiques exacts à l'échelle R (approximés pour affichage)."""
    s = sqrt(3)
    return (
        ("P1 Nord", (0, 0, R)),
        ("P2 Haut-1", (s*R/2, 0, R/2)),
        ("P3 Haut-2", (-s*R/4, 3*R/4, R/2)),
        ("P4 Haut-3", (-s*R/4, -3*R/4, R/2)),
        ("P5 Sud", (0, 0, -R)),
        ("P6 Bas-1", (s*R/2, 0, -R/2)),
        ("P7 Bas-2", (-s*R/2, 0, -R/2)),
    )

def draw(R, filename=None, show=True):
    if R <= 0:
        raise ValueError("Le rayon R doit être strictement positif.")
    import numpy as np
    import matplotlib.pyplot as plt

    pts = coordinates(R)
    coords = [p for _, p in pts]
    fig = plt.figure(figsize=(10.4, 8.4))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor("#f7f9fc")

    # Sphère guide semi-transparente : ne désigne aucune face solide.
    theta = np.linspace(0, 2*pi, 48)
    phi = np.linspace(0, pi, 26)
    x = R * np.outer(np.cos(theta), np.sin(phi))
    y = R * np.outer(np.sin(theta), np.sin(phi))
    z = R * np.outer(np.ones_like(theta), np.cos(phi))
    ax.plot_wireframe(x, y, z, rstride=5, cstride=5,
                      color="#c7d4e0", linewidth=0.40, alpha=0.22)

    for i, j in combinations(range(7), 2):
        a, b = coords[i], coords[j]
        ax.plot([a[0], b[0]], [a[1], b[1]], [a[2], b[2]],
                color="#7b8392", linewidth=1.0, alpha=0.50)

    for idx, (name, point) in enumerate(pts):
        x, y, z = point
        col = "#d28b20" if idx <= 3 else "#1976b7"
        ax.scatter([x], [y], [z], color=col, s=90, edgecolors="white",
                   linewidth=0.8, depthshade=False)
        ax.text(x*1.075 + R*.02, y*1.075 + R*.02, z*1.075,
                name, fontsize=10, color="#17263b")
        ax.plot([0, x], [0, y], [0, z], color=col,
                alpha=0.35, linestyle=":", linewidth=1.1)

    # Le centre O et le barycentre G ne sont pas des sommets.
    ax.scatter([0], [0], [0], color="black", marker="x", s=90)
    ax.text(0, 0, 0, "  O — centre", fontsize=9)
    ax.scatter([0], [0], [R/14], color="darkgreen", marker="D", s=40)
    ax.text(0, 0, R/14, "  G — barycentre", fontsize=9)

    lim = R * 1.20
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.view_init(elev=19, azim=-57)
    ax.set_title("POLYBRUTUS — 7 points sur sphère / 21 distances\n"
                 "Construction exacte : ±60° depuis les pôles",
                 fontsize=13, pad=20)
    fig.tight_layout()
    if filename:
        out = Path(filename)
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=180, bbox_inches="tight")
        print("Image enregistrée :", out.resolve())
    if show:
        plt.show()
    else:
        plt.close(fig)

def main():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--radius", type=float, default=10.0)
    p.add_argument("--save", default=None, metavar="FICHIER.png")
    p.add_argument("--no-show", action="store_true")
    args = p.parse_args()
    draw(args.radius, filename=args.save, show=not args.no_show)

if __name__ == "__main__":
    main()
