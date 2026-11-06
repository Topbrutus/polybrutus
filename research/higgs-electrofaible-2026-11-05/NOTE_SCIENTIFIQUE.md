# Polybrutus × Higgs : holonomie électrofaible, vide neutre et obstruction spectrale

**Rapport de calcul — 2026-11-05.** Auteur du projet : Gabriel St-Pierre. Statut : mathématiques d'un modèle statique 1D, pas théorie des particules prouvée.

## 1. Action et invariance de jauge

Sur sept sites cycliques, soit `phi_i` un doublet de C² de norme 1. Sur l'arête `i→i+1`, soit `S_i∈SU(2)` et `z_i∈U(1)` la phase dans la représentation hypercharge du doublet. L'énergie (a=kappa=lambda=1) est

`E(phi,S,z)=sum_i ||phi_i−z_i S_i phi_(i+1)||² + sum_i (||phi_i||²−1)²`.

Sous les changements locaux `h_i∈SU(2)` et `t_i∈U(1)`, la transformation conjointe

`phi_i→t_i h_i phi_i, S_i→h_i S_i h_(i+1)†, z_i→t_i z_i conj(t_(i+1))`

conserve chaque terme de l'action. Un contrôle négatif laisse les liens fixes après transformation des champs : l'énergie change. La trace du produit des transports `W=(z_0 S_0)...(z_6 S_6)` est invariante par changement de jauge local (conjugaison à la base). Ce résultat algébrique est classique.

## 2. Deux holonomies de référence préenregistrées

Poser `q=(3+4i)/5`, `R=diag(q,conj(q))∈SU(2)`. Sur les six premières arêtes `S=I,z=1`. Sur la dernière :

- **Cas générique SU(2) :** `S=R,z=1`, donc `W=diag(q,conj(q))` et **aucune valeur propre égale à 1**.
- **Cas électromagnétique résiduel :** `S=R,z=q`, donc `W=diag(q²,1)` et `W e_neutre = e_neutre`, avec `e_neutre=(0,1)^T`.

Dans le deuxième cas la configuration `phi_i=e_neutre` à tous les sites a exactement `E=0`, bien que les produits de Wilson SU(2) et U(1) soient **séparément non triviaux**. La neutralité est celle de la composante inférieure du doublet du modèle standard, sous `Q=T3+Y` avec `Y=1/2`. En revanche, le doublet constant chargé `(1,0)` a exactement `E=64/25`.

## 3. Obstruction et minimum exact pour la connexion SU(2) seule

Écrire `theta=arg(q)=atan2(4,3)∈(0,π)`. En diagonalisant l'holonomie, les valeurs propres du laplacien cyclique tordu de 7 arêtes sont

`lambda_(m,±)=2−2 cos[(2πm±theta)/7]`, `m=0,...,6`.

Le plus petit mode donne

`lambda_min=2−2 cos(theta/7)`, et puisque `sum_i||phi_i||²=7`, l'énergie est bornée par `E ≥ 7 lambda_min`. Elle **atteint** cette borne pour les sept doublets neutres `phi_i=exp(−i theta i/7) e_neutre` ; chacun conserve la norme 1. Par conséquent :

**`E_min = 14[1−cos(atan2(4,3)/7)] ≈ 0,122659956426`.**

Une borne inférieure indépendante et élémentaire suit de Cauchy–Schwarz. Avec `Delta_i=phi_i−phi_(i+1)` sur les six premières arêtes et `Delta_6=phi_6−R phi_0`, on a `sum_i Delta_i=(I−R)phi_0`. Or `(I−R)†(I−R)=(4/5)I` et `||phi_0||=1`, donc :

**`E ≥ (1/7)||sum Delta_i||² = 4/35 > 0`.**

Cette preuve d'obstruction est indépendante de la géométrie Polybrutus. Les 14 valeurs propres de l'opérateur discrétisé ont été comparées numériquement à la formule spectrale et une solution explicite vérifie l'égalité.

## 4. Quatre témoins : aucune sélection intrinsèque de Polybrutus

Pour l'attribution initiale (qui reste arbitraire), les énergies dans les quatre fonds prédéfinis sont :

| Géométrie | Identité | SU(2) seule | Compensation électromagnétique | U(1) seule |
|---|---:|---:|---:|---:|
| Polybrutus 3D dans S³ | 15,5 | 15,1 | 15,5 | 15,1 |
| Témoin 4D A | 18 | 17,2 | 17,1466666667 | 17,2 |
| Témoin 4D C | 16 | 17,3636581165 | 17,4363897398 | 17,3636581165 |
| Cube + 1 | 14 | 14,4 | 15,28 | 14,4 |

Toutes ces énergies sont dans des unités du modèle. **Il n'existe aucune relation dérivée ici entre les positions physiques des sites et les points internes Polybrutus.** Réattribuer les points aux sites change l'énergie. Pour le fond identité, la moyenne sur 5040 permutations vaut `65/4=16,25` pour Polybrutus 3D et ses deux témoins 4D partageant les distances, contre 11 pour Cube+1.

Enfin, pour tout jeu arbitraire de sept doublets unitaires `phi_i=(a_i,b_i)`, construire

`F_i=[[conj(b_i),a_i],[-conj(a_i),b_i]]∈SU(2)`, `F_i e_neutre=phi_i`.

Le choix **ajusté a posteriori** `S_ij=F_i F_j†`, `z_ij=1` possède une holonomie totale identité et transporte exactement `phi_j` sur `phi_i`. Il rend `E=0` pour les quatre jeux, y compris Cube+1. Par changement de jauge `h_i=F_i†`, tout couple `(phi,S)` ainsi construit devient le vide constant avec les liens identité.

**Conséquence :** sans mécanisme physique indépendant qui fixe ou contraigne les connexions, l'orientation des sept doublets ne produit pas une signature intrinsèque du Higgs. On ne peut pas « choisir un champ de jauge » après avoir vu Polybrutus, puis invoquer la valeur obtenue comme prédiction.

## 5. Conclusion et prochaine expérience

- **CALCUL :** invariance SU(2)×U(1) conjointe vérifiée dans quatre fonds et quatre configurations; relation exacte entre holonomie et zéro mode; obstruction spectrale exacte dans l'autre fond; contre-exemple de non-identifiabilité.
- **CONNU :** brisure électrofaible et U(1) électromagnétique résiduel sont des propriétés établies du modèle standard. Aucun nouveau mécanisme n'est revendiqué.
- **LIMITES :** 1D statique; pas de 3+1D, pas de terme cinétique temporel, pas de champs dynamiques SU(2) et U(1), pas de couplages physiques, pas de masse du boson, pas de données expérimentales, pas de revue externe, aucune antériorité scientifique établie de notre éventuel exemple particulier.
- **CRITÈRE DE SUITE :** dériver de contraintes indépendantes un placement de sites et une connexion (ou leur action dynamique), puis préenregistrer une observable distinguant Polybrutus des témoins, avec condition d'abandon si cette distinction disparaît sous jauge, permutation ou variation admissible des paramètres.

Ce laboratoire apporte une **limite mathématique falsifiable** à notre hypothèse de travail, non une preuve du boson de Higgs.
