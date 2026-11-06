# Polybrutus — certificat exact d'une réalisation de rang affine 5

**Recherche du 5 novembre 2026 — projet de Gabriel St-Pierre (Topbrutus).** Calcul algébrique et recherche restreinte vérifiés; antériorité mondiale et revue par les pairs **non établies**. Aucun lien physique au boson de Higgs démontré.

## Résultat mathématique

Les sept points originaux de Polybrutus ont les 21 distances **au carré** non étiquetées de multiensemble

`{1×6, (7/4)×2, (9/4)×3, 3×6, (13/4)×2, 4×2}`.

Voici une nouvelle réalisation avec les **mêmes** 21 distances mais un **rang affine exactement 5**. Poser K_ij = 4 d_ij² :

```
0   4   4  12  12   4   4
4   0   4  12  13   9  12
4   4   0  16   9   9  12
12 12  16   0   7   7  12
12 13   9   7   0  13  16
4   9   9   7  13   0   4
4  12  12  12  16   4   0
```

Le calcul exact donne B=8G où G est la Gram des six vecteurs relatifs au premier point :

```
8   4   4   3  -1  -4
4   8   0   7  -1  -4
4   0  24  17   9   4
3   7  17  24   3   0
-1 -1   9   3   8   4
-4 -4   4   0   4   8
```

La preuve n'utilise **aucun arrondi**. Le vecteur entier `w=(4,-5,-5,4,4,0)` est dans le noyau (`Bw=0`), donc rang B ≤5. Le sous-bloc principal aux indices **0,1,2,3,5** possède les mineurs dominants `8,48,1024,3840,14336`, tous positifs. Par le critère de Sylvester, ce sous-bloc est défini positif et de rang 5. Comme le noyau a une composante non nulle sur l'indice omis 4, B est **positive semi-définie de rang 5**. D'où l'existence exacte de sept points en R5 réalisant ces 21 distances. Contrôle : les six mineurs principaux de taille 5 de B sont `0,14336,14336,22400,22400,14336` et `det B=0`.

**Cosphericité :** ce témoin de rang 5 est sur une sphère de rayon au carré **43/28**, et non de rayon 1. La condition de compatibilité provient de `w·diag(B)=0` ; la valeur du rayon se calcule exactement avec le sous-bloc défini positif.

## Méthode de découverte et limitations

Le témoin 4D précédent est conservé pour ses six premiers sommets (15 distances fixes). Les 6 distances du septième sommet sont réattribuées par toutes les permutations distinctes de `[16,12,12,4,4,12]`, soit **60 affectations**. La recherche rationnelle exhaustive **dans cette famille restreinte** identifie deux candidats PSD : l'original de rang 4, et le nouveau rang 5. Ce dénombrement ne donne **aucune** unicité globale parmi toutes les affectations des 21 distances.

Les réalisations de rang 3, 4 et 6 sont déjà documentées. Nous possédons désormais une réalisation **3, 4, 5 et 6** pour le même multiensemble. Ces résultats n'impliquent pas un nouveau Higgs ni un champ à cinq composantes.

## Reproduction

- [Vérificateur autonome](preuve_minimale_5d.py), ne requérant que Python 3 et sa bibliothèque standard.
- Un laboratoire complet construit dans la conversation contient le programme détaillé, 32 tests unitaires, un audit SymPy indépendant, des données JSON et 9 empreintes SHA-256 contrôlées après réextraction de l'archive.
- La **note GitHub et le vérificateur minimal** sont enregistrés ici; cette opération ne prétend pas déposer l'archive complète sur GitHub.

**Prochaine étape :** comparer ce cas de reconstruction des distances non étiquetées aux travaux de géométrie des distances, solliciter une relecture de la preuve et garder le modèle électrofaible comme programme physique séparé.
