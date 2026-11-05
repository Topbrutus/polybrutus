# NOTE SCIENTIFIQUE — Polybrutus, sphère du Higgs et transport de jauge

**2026-11-05 — Gabriel St-Pierre (Topbrutus), assistance computationnelle**  
**État :** résultats de mathématiques et de covariance de jauge vérifiés dans un modèle discret; **aucune identification à une particule**. La nouveauté dans la littérature n'a pas été établie.

## 1. L'hypothèse précise

Une orientation d'un doublet complexe de Higgs possède quatre coordonnées réelles. À norme non nulle fixée, les orientations sont paramétrées par `S³ ⊂ R⁴`. Le projet Polybrutus donne sept points exacts sur une sphère, en 3D ou en 4D. L'analogie « sept états sur une sphère » est mathématiquement possible mais ne prouve ni que le champ de Higgs comporte sept constituants ni que les points constituent des entités physiques.

Le vide du champ de Higgs possède une orbite de jauge; une transformation de `SU(2)` agit transitivement sur les doublets unitaires. Pour un champ de Higgs local, le nombre de degrés de liberté physiques n'est pas le nombre de directions de cette orbite. Dans le modèle électrofaible, trois modes angulaires sont absorbés par les bosons de jauge (mécanisme de Higgs), et un mode radial demeure.

## 2. Lemme SU(2) exact

Écrire un doublet normalisé `ψ=(a,b) ∈ C²`, avec `|a|²+|b|²=1`. Définir :

```text
G(ψ) = [ a       −conj(b) ]
       [ b        conj(a) ]
```

Alors `G†G=I`, `det G=1` et `G(ψ)(1,0)^T=ψ`. Le choix local `h_i=G(ψ_i)†` donne **ψ'_i=(1,0)^T** pour chacun des sept sites. Ce fait ne représente pas une nouvelle découverte de physique : il traduit la transitivité de `SU(2)` sur `S³`.

Sans liaisons de jauge, les différences `||ψ_i−ψ_j||²` varient sous les transformations indépendantes des sites. Si on les aligne toutes, ces différences deviennent nulles : **les 21 distances ne sont pas intrinsèquement des observables électrofaibles**.

Précision : des transformations indépendantes sont légitimes pour des **sites spatiaux distincts**; les sept sites n'ont pas encore de localisation physique définie. Une transformation commune à tous les points préserverait au contraire leurs distances.

## 3. L'invariant correctement transporté

Sur un graphe complet de sept sites, introduire `U_ij ∈ SU(2)` et la loi :

```text
ψ_i'  = h_i ψ_i
U_ij' = h_i U_ij h_j†
I_ij  = ||ψ_i − U_ij ψ_j||²
```

L'identité `ψ_i'−U_ij'ψ_j'=h_i(ψ_i−U_ijψ_j)` prouve rigoureusement `I'_ij=I_ij`. Les liaisons `U_ij` supposent ici un transport parallèle discret non dérivé d'un champ électrofaible de l'espace-temps.

Une grandeur à trois sites peut être formée :

```text
T_U = Σ_{i<j<k} I_ij I_ik I_jk.
```

C'est un **descripteur discret invariant de jauge**, **pas** une interaction physique fondamentale ni un terme spontané du lagrangien du modèle standard. Son caractère invariant provient du transport; sa signification physique demeure entièrement ouverte.

Le barycentre transporté vers un site de référence 0,

```text
B_0 = (1/7) Σ_i U_0i ψ_i,  avec U_00=I,
```

se transforme comme `B_0' = h_0 B_0`. Sa norme au carré est donc invariante de jauge **une fois la connexion spécifiée**. Le barycentre ordinaire `(1/7)Σψ_i` ne l'est pas sous des transformations locales indépendantes.

## 4. Résultats exacts

Normes des doublets unitaires, matrice de Gram, cosphéricité, matrices `SU(2)`, sommes sur les 21 liens et 35 triangles ont été contrôlées. Tous les quatre jeux testés sont cosphériques et normalisés sur `S³`. Le témoin `4D_B` (non cosphérique) est volontairement exclu du test de vide à norme fixe.

| Jeu | `ΣI_ij` pour `U=I` | `T_U` pour `U=I` | `||B||²` pour `U=I` | Après changement de jauge cohérent | Connexion pure `U_ij=G_i G_j†` |
|---|---:|---:|---:|---|---|
| Polybrutus 3D plongé dans 4D | 195/4 | 24435/64 | 1/196 | valeurs **inchangées** | `ΣI=0`, `T=0`, `||B||²=1` |
| 4D A | 195/4 | 23835/64 | 1/196 | valeurs **inchangées** | `ΣI=0`, `T=0`, `||B||²=1` |
| 4D C | 195/4 | 24273/64 | 1/196 | valeurs **inchangées** | `ΣI=0`, `T=0`, `||B||²=1` |
| Cube six faces + septième point | 33 | 135 | 16/49 | valeurs **inchangées** | `ΣI=0`, `T=0`, `||B||²=1` |

**Une nuance cruciale :** `U=I` et `U_ij=G_i G_j†` sont *deux choix différents de connexion* si on conserve les mêmes `ψ_i`; ils ne sont pas simplement deux représentations de jauge d'une même paire (champ, connexion). Leur différence dans `ΣI` ne contredit pas la covariance. Les deux connexions sont plates (`tr(U_ij U_jk U_ki)=2` sur les 35 triangles). La platitude seule ne fixe pas la variation covariante du champ.

**Test supplémentaire :** en appliquant réellement `ψ'_i=h_i ψ_i` et `U'_ij=h_i U_ij h_j†`, le barycentre nu passe de `||\bar ψ||²=1/196` à `1` pour les trois jeux Polybrutus, et de `16/49` à `1` pour Cube+1. Cependant, `||B_0||²` reste inchangé. Pour les jeux 3D et Cube+1, toutes les vérifications ont été développées en matrices symboliques; pour A et C, le résultat est certifié par le lemme SU(2) avec normes, matrices de distance et frames exacts, et par au moins un contrôle direct d'arête.

## 5. Ce que le test n'établit pas

- Il ne construit **pas** un modèle du boson de Higgs, un état lié, une prédiction de masse, ni une nouvelle interaction.
- Il ne fournit **pas** les positions d'espace-temps des sept sites, les valeurs de `g` et `g'`, la connexion électrofaible tirée des équations du mouvement ou un réseau physique validé.
- Une somme de triangles invariante ne devient pas spontanément un terme de lagrangien : il faudrait justifier l'opérateur, sa dimension, ses coefficients, son comportement sous renormalisation et son pouvoir prédictif.
- Le champ de Higgs est une **section de fibré associée**, pas une boule physique; `S³` est une description d'orbite de valeurs à norme constante.
- Le rapport `4/7` du modèle Cube+1 provient du choix angulaire 60°, pas d'une loi électrofaible. Le fait que la *norme du barycentre nu* change sous jauge interdit de l'identifier comme grandeur physique sans transport parallèle.
- Les résultats du témoin B et les preuves 3D–4D antérieures ne sont ni altérés ni remplacés.

## 6. Critère prospectif pour la prochaine expérience

**But :** déterminer si une géométrie à sept sites améliore un modèle physique de référence, dans une théorie dont la connexion et les couplages sont fixés indépendamment de la cible recherchée.

1. Définir les sites `x_i` et les chemins ou le maillage avant de regarder les résultats.
2. Choisir un modèle de connexion `SU(2)_L × U(1)_Y` dont l'origine physique est justifiée (et sa règle de transport).
3. Construire une discrétisation contrôlée du terme cinétique `|D_μΦ|²`, avec unités et facteur d'échelle.
4. Comparer les mêmes observables invariantes avec des jeux témoins et un modèle sans structure Polybrutus.
5. Préenregistrer un test pouvant **réfuter** le modèle et rendre publics les échecs.

**Condition d'abandon de l'interprétation physique :** si aucun opérateur ou paramètre de jauge physiquement motivé ne distingue notre construction des témoins, conserver Polybrutus comme objet mathématique, sans le présenter comme une description de particule.

## 7. Reproduction et provenance

Dans le dossier `research/higgs-jauge-2026-11-05/` du dépôt :

```bash
python -m pip install sympy==1.14.0
python -B higgs_gauge_links.py
python -B -m unittest -v test_higgs_gauge_links.py
```

**17 tests passent** en local dans Python, SymPy 1.14.0; le journal détaillé, JSON de résultats et manifeste SHA-256 accompagnent cette note. Source originale des distances : `../higgs-variationnel-2026-11-05/polybrutus_higgs_lab.py` et les documents de référence de Polybrutus sur `main`.

Références de fond : formulation usuelle du doublet de Higgs et du mécanisme de Higgs (cours CERN : https://indico.cern.ch/event/763013/contributions/3358697/attachments/1813182/2962418/Higgs-1.pdf). La règle de transformation des liaisons est celle du transport parallèle en théorie de jauge sur un réseau (Wilson links). Cette note ne revendique aucune nouveauté sur ces outils théoriques établis.
