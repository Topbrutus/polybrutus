# Polybrutus — reconstruction par 21 distances non étiquetées (3D / 4D)

**Étude complémentaire à [Polybrutus v0.1.0](../../README.md).** Auteur du projet : Gabriel St-Pierre (Topbrutus). **Statut : recherche computationnelle reproductible, non évaluée par les pairs.** Aucune revendication de nouveauté mondiale, d'optimalité ou de résultat physique.

La configuration initiale de sept points possède 21 distances par paires. Avec un rayon de référence de 1, leurs **carrés** forment le multiensemble :
`{1 (×6), 7/4 (×2), 9/4 (×3), 3 (×6), 13/4 (×2), 4 (×2)}`.

Cette étude documente deux assertions **de portée limitée** :
1. Deux énumérateurs exacts de réalisations dans **R³**, sans supposer d'avance une sphère commune, produisent 27 complétions étiquetées, toutes congruentes à la configuration originale; aucune réalisation plane ou alignée n'a été trouvée.
2. Un **témoin exact en R⁴**, de rang affine 4, possède le même multiensemble de distances mais n'est pas congruent à Polybrutus (rang affine 3). Sa matrice de Gram contient un mineur 4×4 de déterminant `15/16`.

**Ce qui reste ouvert :** revue spécialisée de la complétude des énumérateurs, recherche d'antériorité approfondie, validation externe et preuve analytique éventuelle. Le simple fait d'exécuter ces programmes ne remplace pas une relecture mathématique.

## Reproduction

Testé le 5 novembre 2026 sous **Python 3.12.10** et **SymPy 1.14.0** sur Windows. Dans ce dossier :

```bash
python -m pip install sympy==1.14.0
python -B AUDIT_INDEPENDANT_COORDONNEES_2026_11_05.py
python -B polybrutus_contre_exemple_4d_exact_2026_11_05.py
python -B polybrutus_unicite_sans_sphere_2026_11_05.py
python -B polybrutus_independant_cayley_menger_2026_11_05.py
python -B polybrutus_controle_sympy_global_2026_11_05.py
```

Les deux énumérateurs R³ utilisent uniquement la bibliothèque standard; les trois audits complémentaires utilisent SymPy. `laboratoire_spectral.py` est ici une **dépendance géométrique minimale**, conservant les coordonnées exactes et `d2e` du laboratoire d'origine, sans ses données expérimentales annexes.

- [Note scientifique / Scientific review note](NOTE_RELECTURE_SCIENTIFIQUE_2026_11_05.md)
- [Contrôles d'exécution](VERIFICATION_2026_11_05.md)

**Zenodo :** ce répertoire ne constitue pas une nouvelle publication Zenodo et ne demande aucun nouveau DOI. Toute future publication complémentaire devra être explicitement validée.
