# Polybrutus : neuf reperes et realisations en dimensions 3, 4 et 6

**2026-11-05. Projet de Gabriel St-Pierre (Topbrutus). Recherche en geometrie euclidienne exacte. Non soumise a revue par les pairs; nouveaute mondiale non etablie. Aucun lien physique avec le boson de Higgs demontre.**

## Motivation: deux observations d'un dessin
Le croquis fourni dans la conversation montre une figure numerotee de 1 a 9 ainsi qu'une transition symbolique 3D, 4D et 5D. Cette image n'est pas une demonstration. L'etude ci-dessous examine une construction mathematique precise, sans identifier le croquis a Polybrutus.

## Neuf points de reference (7 sommets + O + G)
Les sept sommets sont ceux de SPECIFICATION.md (R=1); O=(0,0,0) est le centre de la sphere et G=(0,0,1/14) leur barycentre. O et G ne sont pas de nouveaux sommets.
Les neuf references sont distinctes et forment 36 paires de distances. Leur histogramme exact des distances au carre est :

| d² | Multiplicite |
|---|---:|
| 1/196 | 1 |
| 169/196 | 1 |
| 183/196 | 3 |
| 1 | 13 |
| 211/196 | 2 |
| 225/196 | 1 |
| 7/4 | 2 |
| 9/4 | 3 |
| 3 | 6 |
| 13/4 | 2 |
| 4 | 2 |

Somme des 36 distances² : **12293/196**. Rang affine des 9 references : 3. La distance OG²=1/196 n'est pas independante de la cosphericite unitaire et du multiensemble initial, car sum_{i<j}d_ij²=49(1-|G-O|²)=195/4.

## Realisation nouvelle en dimension 6: certificat exact
Les sept distances originales ont le multiensemble d² = {1×6, 7/4×2, 9/4×3, 3×6, 13/4×2, 4×2}.
Pour K_ij = 4d_ij², une realisation euclidienne est certifiee par la matrice 6x6 (8G)_ij = K_{1,i+1}+K_{1,j+1}-K_{i+1,j+1}. Si 8G est positive definie, la realisation possede le rang affine 6.

```text
 0  4 12 12 16  9 12
 4  0  7 12 13  4 13
12  7  0  9  4 12  7
12 12  9  0  4  9  4
16 13  4  4  0 12  4
 9  4 12  9 12  0 16
12 13  7  4  4 16  0
```

Cette matrice K est symetrique, de diagonale nulle, et possede exactement les multiplicites requises pour les 21 distances. **Les six mineurs principaux dominants de 8G sont** : `8, 111, 1560, 3489, 6324, 4688`, tous strictement positifs. Le critere de Sylvester prouve donc l'existence exacte d'une realisation en R^6 de rang affine 6. Cette configuration n'est congruente ni au temoin de rang 3, ni au temoin de rang 4. L'antériorité dans la litterature n'a pas ete etablie.

## Quelle interpretation du croquis 3D + 3D -> 4D/5D ?
Pour deux vrais sous-espaces vectoriels U et V de dimension 3, dim(U+V)=6-dim(U inter V). On obtient 4 dimensions si l'intersection est de dimension 2, 5 si elle est de dimension 1, 6 si elle est de dimension 0. Cette identite connue n'interprete pas physiquement le croquis.

## STATUT D5: question OUVERTE
Aucune matrice de Gram PSD de rang **5** ayant ce multiensemble de 21 distances n'est actuellement certifiee ni interdite. Les certificats de rang 3, 4 et 6 **n'impliquent pas** un certificat de rang 5. Il faut trouver une affectation des 21 longueurs aux 21 paires telle que la Gram soit PSD et de rang exactement 5, ou demontrer que c'est impossible.

## Preuves de calcul et limites
Un laboratoire independant en Python standard (arithmetique Fraction, sans flottants) a verifie les 3 rangs PSD par tous les mineurs principaux, les 36 distances des 9 points, les 6 mineurs de Sylvester de rang 6, et 20 tests unitaires. Archive locale ZIP avec code, journal, JSON et 6 SHA-256 verifiés et re-extraits dans un dossier vierge. **Ce fichier GitHub est une note compacte; les sources de l'archive ne sont pas pretendues publiees ici.**

Le doublet de Higgs physique possede 4 composantes reelles internes; les rangs affines 5/6 de geometries euclidiennes ne prouvent aucune nouvelle particule ni extension du Modele standard. Prochaine action : resoudre la question exacte D5 et rechercher l'anteriorite des certificats.
