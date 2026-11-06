# Note de recherche — Polybrutus sur la sphère de vide du Higgs

**Date :** 2026-11-05  
**Auteur du projet géométrique :** Gabriel St-Pierre (Topbrutus)  
**Statut :** exploration analytique et computationnelle; **non soumise à évaluation par les pairs**.  
**Objet :** tester l'intuition de sept directions géométriques encadrées par une sphère liée à la structure du champ de Higgs, sans extrapoler à une masse de particule.

## 1. Hypothèse testée et distinction des objets

**Source géométrique :** Polybrutus définit sept points, 21 distances non étiquetées avec six valeurs au carré `M={1×6, (7/4)×2, (9/4)×3, 3×6, (13/4)×2, 4×2}`. Les données comprennent Polybrutus 3D (rang affine 3), le témoin 4D A déjà connu (rang 4), et deux témoins 4D B/C issus d'un balayage exploratoire local, certifiés ensuite exactement.

**Hypothèse :** comparer les configurations cosphériques au sous-ensemble de valeurs de vide d'un doublet complexe de Higgs. C'est une **mise en correspondance mathématique**, pas une déduction du modèle standard. L'orientation d'un doublet de Higgs en chaque point de l'espace-temps n'est pas une nouvelle particule.

**Contrôles :** modèle cube `6+1` de la note antérieure, et 4D B **non cosphérique** comme contre-exemple à l'implication « distances M => sphère ».

## 2. Potentiel variationnel exact du Higgs

Pour un doublet complexe `Φ` contenant quatre composantes réelles, écrire la norme `Φ†Φ=(v+h)^2/2` lors d'une variation radiale et choisir, par convention, le potentiel

```
V(Φ) = λ [Φ†Φ − v²/2]²   (λ > 0, v > 0).
```

Si `u∈S³⊂R⁴` est une orientation normalisée (après identification réelle/complexe), `Φ=(v+h)u/√2`. Il en découle sans approximation :

```
V(h) = (λ/4)(2 v h + h²)²
     = λ v² h² + λ v h³ + (λ/4) h⁴.
V(0)=0,  dV/dh|0=0,  d²V/dh²|0=2 λ v².
```

Pour le champ radial canoniquement normalisé, la relation habituelle au niveau arbre est `m_H²=2 λ v²` en unités naturelles. La masse **n'est pas prédite par Polybrutus** : il faut connaître `λ` et `v` de façon indépendante. Les **trois directions tangentielles** au lieu de vide sont plates pour ce potentiel; elles correspondent aux degrés de liberté Goldstone, qui sont absorbés dans les bosons de jauge de la théorie électrofaible.

**No-go élémentaire :** toutes les orientations de norme constante ont la même valeur de `V`. Le potentiel standard isolé ne sélectionne ni la configuration Polybrutus ni le nombre sept.

## 3. Théorème simple — indistinguabilité par énergie additive à deux corps

Si deux réalisations de sept points ont le même *multiensemble* de 21 distances au carré `M`, alors pour toute fonction `f` définie sur `M` :

```
E_f(P) = Σ_{i<j} f(d_ij²)
       = Σ_{x∈M} f(x)
       = E_f(Q).
```

La preuve est une réindexation de la somme, sans hypothèse d'espace ambiant, de sphère ou de physique.

Pour `f(x)=1/x`, les quatre configurations donnent exactement :

```
E2 = 6 + 8/7 + 4/3 + 2 + 8/13 + 1/2 = 6329/546.
```

**Conséquence :** une fonctionnelle intégralement additive par paires, ne dépendant que des distances, ne peut pas privilégier Polybrutus parmi ces réalisations. Une énergie d'interaction physique spécifique ne découle pas de cette égalité.

## 4. Un descripteur géométrique à trois corps distingue les classes

Définir une quantité *exploratoire* et symétrique sous les permutations de sommets :

```
T(P)=Σ_{i<j<k} d_ij² · d_ik² · d_jk².
```

C'est un descripteur géométrique, **pas une nouvelle force ni une nouvelle énergie fondamentale**. Contrairement à `E_f`, il utilise la manière dont les distances se regroupent dans les triangles.

| Ensemble | Rang affine | Sur sphère | E2 | T | tr(Gram centré³) |
|---|---:|---|---|---|---|
| Polybrutus 3D | 3 | oui, R=1 | 6329/546 | 24435/64 | 4629285/87808 |
| 4D A | 4 | oui, R=1 | 6329/546 | 23835/64 | 5234925/87808 |
| 4D B | 4 | **non** | 6329/546 | 26817/64 | 3972321/87808 |
| 4D C | 4 | oui, R=1 | 6329/546 | 24273/64 | 4721139/87808 |

Le troisième moment spectral est `tr(G_c³)` pour `G_c=-1/2 H D² H` et `H=I−11ᵀ/7`. Tous les éléments affichés sont des fractions exactes calculées avec Python `fractions.Fraction`.

**Limite :** le choix de `T` est exploratoire. Pour lui attribuer un rôle physique, il faudrait dériver ce terme d'un modèle invariant de jauge, assorti d'un couplage dimensionné et de prédictions indépendantes. Il ne suffit pas de constater quatre valeurs différentes.

## 5. Condition sphérique : trois admissibles, une exclue

Pour chaque jeu, le laboratoire reconstruit la matrice de Gram par les distances :

`G_ij=(K_0i+K_0j−K_ij)/8`, avec `K=4d²`.

Il prouve que les matrices sont PSD de rang 3 ou 4 par mineurs principaux positifs et complément de Schur nul. Il calcule la sphère circonscrite des `d+1` premières positions affinement indépendantes, puis vérifie les autres par résidus exacts.

- Polybrutus 3D : rayon au carré 1. Embeddable dans `R⁴` en ajoutant une quatrième coordonnée 0; la configuration reste de rang 3, donc appartient à un sous-ensemble `S²` de `S³`.
- A et C : rayon au carré 1 et rang affine 4; normalisées autour de leurs propres centres, elles donnent des septuplets sur `S³`.
- B : première sphère circonscrite de rayon au carré `1613/1280`, mais résidus non nuls pour d'autres sommets : **aucune sphère commune**. Ne pas l'utiliser comme sept états de norme identique sans changer les distances.

L'appartenance à `S³` n'implique pas que les distances entre sept valeurs hypothétiques du champ soient observables physiques.

## 6. Contre-test variationnel : répulsion inverse-carrée

Choix illustratif, sans paramètres calibrés : `E2=Σ 1/d_ij²` sur sept vecteurs unitaires de `S³`. La valeur d'énergie est identique, mais l'équilibre contraint exige que la projection tangente du gradient en chaque point soit nulle. Calcul exact :

| Configuration sphérique | Somme des normes² des gradients tangents |
|---|---:|
| Polybrutus 3D | `34291555126/5554571841` |
| A | `12874444982806/449920319121` |
| C | `159152051260/21424777101` |

Toutes les valeurs sont **strictement positives**, donc aucun de ces trois échantillons n'est un point stationnaire de cette répulsion simple, lorsque tous les points peuvent varier indépendamment sur la sphère.

**Attention :** ce test ne falsifie pas le modèle standard du Higgs : la répulsion entre sept points est un choix exploratoire, absent du potentiel standard. Il falsifie seulement l'affirmation précise « l'énergie E2 choisie est stationnaire sur Polybrutus ».

## 7. Contrôle « six faces du cube + référence »

En prenant `P0=e0` et `P_(i,±)=cos(θ)e0±sin(θ)e_i` pour `i=1,2,3`, on a sept points unitaires. Le barycentre vaut :

`Gθ=((1+6cosθ)/7)e0`.

À `θ=60°`, `G=(4/7)e0` exactement. Mais le multiensemble des distances au carré est `{1×6, (3/2)×12, 3×3}`, **trois classes**, non les six familles Polybrutus. À `θ=90°`, le barycentre devient `e0/7`. Le **4/7 n'est pas une prédiction du Higgs ni une obligation géométrique**; il découle du choix `θ=60°`.

## 8. Point bloquant : symétrie de jauge et observable

La variété des orientations unitaires du doublet de Higgs possède trois paramètres, mais les orientations du vide sont reliées par `SU(2)_L×U(1)_Y` dans le modèle standard. On ne peut pas promouvoir les distances entre composantes de sept vecteurs au rang d'observables physiques simplement parce qu'elles sont invariantes sous rotations globales `O(4)`.

Il manque dans l'hypothèse : **les positions d'espace-temps `x_i`**, la **connexion de jauge `D_μ`**, un **transport parallèle** permettant de comparer des orientations en différents points, et une **quantité invariante de jauge** (ou un observable physique directement calculable). Sans cela, les sept points sont un choix de coordonnées/échantillonnage et non une configuration sélectionnée par le modèle électrofaible.

L'image des dérivées tensorielles est pertinente à ce stade; Einstein–Hilbert, Riemann et surfaces minimales ne sont pas intégrés sans mécanisme explicite. Il serait incorrect de coller leurs formules au modèle pour améliorer une coïncidence numérique.

## 9. Hypothèse suivante et condition d'abandon

**Candidat de travail, non adopté comme loi physique :** choisir un espace-temps de test, identifier ce que représenteraient les sept emplacements `x_i`, introduire une connexion électrofaible définie indépendamment des résultats et calculer une quantité gauge-invariante (par exemple un corrélateur avec transport). Comparer avec un modèle témoin sans structure à sept points; annoncer les paramètres et le critère de comparaison *avant* d'examiner une cible expérimentale.

**Test minimal :** une même configuration soumise à des transformations de jauge indépendantes doit laisser l'observable proposé inchangé. Si l'« invariant » varie, abandonner son interprétation physique ou reformuler une quantité réellement invariante. Si tous les observables candidats restent égaux pour l'ensemble des configurations, alors la géométrie à sept points n'apporte pas de distinction physique dans ce cadre.

**Interdits pour cette étude :** retoucher la V1 figée du boson de Brutus, choisir un coefficient pour reproduire `125 GeV`, attribuer des unités de masse à une longueur sans dérivation, présenter une identité après-coup comme une prédiction.

## 10. Reproductibilité et sources

```bash
python -B polybrutus_higgs_lab.py
python -B -m unittest -v test_polybrutus_higgs_lab.py
```

10 tests verts lors de la passe locale du 5 novembre 2026. Aucune dépendance externe. Les fichiers de preuve et les SHA-256 sont fournis dans ce dossier.

**Sources de géométrie :**
- Polybrutus : https://github.com/Topbrutus/polybrutus/blob/main/SPECIFICATION.md
- Étude de reconstruction : https://github.com/Topbrutus/polybrutus/tree/main/research/distance-reconstruction-3d-4d
- Certificat local des classes B/C : `polybrutus_certificat_trois_classes_4d_2026_11_05.py`, déjà produit dans une session précédente (source archivée séparément).
- Contrôle cube : note et script locaux du 5 novembre 2026.

**Sources de physique :**
- Cours du CERN sur le mécanisme de Higgs (quatre degrés de liberté, trois Goldstone) : https://indico.cern.ch/event/763013/contributions/3358697/attachments/1813182/2962418/Higgs-1.pdf
- Présentation des rapports m_H²=2 λv² (convention de λ explicitée) : https://link.springer.com/article/10.1140/epjc/s10052-014-2744-3

**Conclusion :** géométrie exacte; discrimination à trois corps démontrée; absence de sélection par le potentiel radial du Higgs; **aucun nouveau résultat physique validé**.
