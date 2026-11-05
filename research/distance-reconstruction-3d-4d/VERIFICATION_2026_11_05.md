# Audit local des calculs — 2026-11-05

**Environnement observé :** Windows / Python 3.12.10 / SymPy 1.14.0.

- Manifeste SHA-256 du dossier de relecture d'origine : **14/14 fichiers conformes**.
- `AUDIT_INDEPENDANT_COORDONNEES_2026_11_05.py` : **PASS**, multiensembles des 21 distances identiques, 5 040 permutations sans congruence, rangs 3 et 4, perturbation négative détectée.
- `polybrutus_contre_exemple_4d_exact_2026_11_05.py` : **PASS**, toutes les 21 distances et mineur de Gram exact `15/16`.
- `polybrutus_unicite_sans_sphere_2026_11_05.py` : **PASS**, rang 3 : 27 configurations compatibles / 27 congruentes ; rangs 2 et 1 : aucune solution.
- `polybrutus_independant_cayley_menger_2026_11_05.py` : **PASS**, mêmes dénombrements (27/27/0/0) avec déterminants de Cayley–Menger.
- `polybrutus_controle_sympy_global_2026_11_05.py` : **PASS**, validation algébrique des 27 matrices exactes.

**Portée du test :** les cinq programmes ont été exécutés sur le PC de travail avec le module complet de laboratoire d'origine. Le module présent dans ce sous-dossier est une version réduite ne contenant que `EXACT_POINTS` et `d2e`; ses deux définitions proviennent de la même source et doivent être retestées par tout relecteur sur un clone propre.

**Statut scientifique :** revue externe NON EFFECTUÉE; priorité scientifique INCONNUE. La complétude mathématique des réductions par symétrie reste un élément de relecture, malgré deux implémentations concordantes.

**Non-régression :** aucune modification des fichiers historiques à la racine du dépôt, et aucune opération sur Zenodo.
