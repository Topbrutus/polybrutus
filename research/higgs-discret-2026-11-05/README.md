# POLYBRUTUS × HIGGS — laboratoire discret de jauge à sept sites

**2026-11-05 — Recherche exploratoire, Gabriel St-Pierre (Topbrutus)**  
**Statut :** calcul de modèle-jouet exact ; **aucune prédiction de masse, aucune validation du Higgs, aucune nouvelle particule démontrée**.

## Question falsifiable

Avec des sites d'espace **distincts** des orientations internes du doublet de Higgs, la géométrie Polybrutus peut-elle donner des valeurs différentes à une action discrète invariante de jauge ? Et les différences survivent-elles à un contre-test de réattribution des sites ?

## Choix figés pour la passe

- **Topologie spatiale :** sept sites sur une boucle fermée 1D, voisins (0,1), (1,2), …, (6,0), sites repérés sur un heptagone régulier dans un plan, avec longueur de chaque arête spatiale `a=1`.
- **Valeurs internes :** doublets complexes `phi_i∈C²`, normalisés `phi†phi=1`, obtenus par plongement des configurations exactes Polybrutus 3D / 4D A / 4D C et Cube+1 dans `S³`.
- **Énergie *statique* et non normalisée aux unités physiques :**

  `E[phi,U] = (kappa/a) * Σ_(i,i+1) ||phi_i − U_(i,i+1) phi_(i+1)||² + a * Σ_i lambda (phi_i†phi_i − v²/2)²`.

- **Paramètres de démonstration :** `a=kappa=lambda=g=g'=1`, `v=sqrt(2)`, hypercharge du doublet `Y=1/2`. Ce sont des **conventions sans unités** : aucune n'est une mesure physique, un ajustement aux données ou une prédiction.
- **Fond témoin :** tous les transports `U=I₂`. **Fond à holonomie non triviale :** tous `I₂` sauf l'arête `6→0` où `U=qR`, `q=(3+4i)/5`, `R=[[3/5,-4/5],[4/5,3/5]]∈SU(2)`. L'holonomie est fixée dans le code, et n'est **pas dérivée d'une solution des équations électrofaibles**.
- **Correspondance :** la liste ordonnée des 7 sommets est attribuée dans l'ordre aux 7 sites. Cet appariement **n'est pas imposé par la géométrie** et fait partie des hypothèses.

Les paramètres ont été fixés pour cette expérience **avant l'évaluation des énergies de cette passe** ; ils n'ont pas été préenregistrés auprès d'un tiers avant les travaux antérieurs et n'ont aucun statut d'étalonnage expérimental.

## Reproduire

Dans le dépôt GitHub, disposer des répertoires frères `../higgs-variationnel-2026-11-05` et `../higgs-jauge-2026-11-05` non modifiés. Dans l'archive ZIP, les deux fichiers Python sources sont également inclus en tant que dépendances figées.

```bash
python -m pip install sympy==1.14.0
cd higgs-discret-2026-11-05
python -B higgs_lattice_7.py
python -B -m unittest -v test_higgs_lattice_7.py
sha256sum -c MANIFEST_SHA256.txt
```

Le programme émet un JSON de résultats déterministes, et les tests vérifient les invariants exacts, les deux connexions, la covariance locale `SU(2)×U(1)` et les contrôles négatifs. Le fichier `TESTS_2026_11_05.log` est la trace de la passe à l'origine du manifeste, avec sa durée mesurée (pas une durée reproductible).

## Résultats et limites

Avec les transports identités (énergie potentielle nulle au vide), les actions exactes valent `31/2`, `18`, `16`, `14` pour Polybrutus 3D, A, C et Cube+1. Le contre-témoin non cosphérique B est exclu des doubles unités. La connexion non triviale modifie ces valeurs et produit une holonomie `tr W = 18/25 + 24i/25`, invariant sous jauge, distinct du témoin `tr W = 2`.

**Contrôle négatif crucial :** permuter les champs attribués aux sites 3 et 4 change l'énergie de Polybrutus 3D de `31/2` à `37/2` sans changer le multiensemble des 21 distances. Il est donc impossible de tirer une grandeur physique univoque de ces résultats sans justifier indépendamment la topologie, les positions, l'appariement et la connexion.

Un anneau spatial statique à sept sites n'est **pas** le lagrangien électrofaible dynamique en `3+1` dimensions : il ne contient ni équations de mouvement temporelles, ni dynamique de jauge, ni plaquettes locales, ni mécanisme prédisant `m_H`. Le calcul de la courbure du potentiel `m_toy²=4` est une identité du choix `v,lambda`, **pas une masse de particule**.

Voir `NOTE_SCIENTIFIQUE.md` et `PROVENANCE.md` pour les hypothèses, les limites, les tests et l'étape suivante.
