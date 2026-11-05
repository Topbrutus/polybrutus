# POLYBRUTUS–HIGGS — laboratoire variationnel exploratoire

**5 novembre 2026 — Topbrutus / Gabriel St-Pierre**  
**Statut :** calculs géométriques exacts vérifiés, **aucune identification physique établie**.  
**Portée :** quatre jeux de 7 points (Polybrutus 3D, témoins 4D A/B/C) + contrôle « six faces du cube + référence » sur S³.  
**Règle :** aucune modification de la formule du boson de Brutus V1 figée; aucune nouvelle masse prédite.

## Exécuter

Python 3.10+; **aucune bibliothèque externe** nécessaire.

```bash
python -B polybrutus_higgs_lab.py
python -B -m unittest -v test_polybrutus_higgs_lab.py
```

Les résultats archivés de la passe locale sont dans `RESULTATS_2026_11_05.txt` et `TESTS_2026_11_05.log`. Le fichier `MANIFEST_SHA256.txt` contrôle les fichiers de travail sans se compter lui-même.

## Conclusions vérifiées

- Pour toutes les quatre configurations Polybrutus, les **21 distances au carré** ont le même multiensemble, même lorsque B n'est pas cosphérique.
- Toutes les énergies **purement additives par paires** `sum(f(d_ij²))` coïncident donc, quelle que soit la fonction `f` définie sur ces distances. Exemple : `sum(1/d_ij²)=6329/546`.
- Le **descripteur à trois corps** `T=sum_{i<j<k} d_ij²*d_ik²*d_jk²` varie : 3D `24435/64`, A `23835/64`, B `26817/64`, C `24273/64`.
- La classe 4D B n'a pas de sphère commune et ne peut **pas** être directement placée sur la sphère de vide S³ sans transformation qui changerait les distances.
- Le potentiel de Higgs standard, à lui seul, ne sélectionne aucune des configurations sphériques : il ne dépend que de la norme du champ. Autour du vide, il y a une direction radiale massive et trois directions tangentielles plates pour le potentiel.
- Aucun des trois jeux Polybrutus sphériques n'est stationnaire pour **l'énergie géométrique illustrative** `sum(1/d²)` sur sa sphère de rayon 1. Celle-ci **n'est pas** une interaction physique dérivée du modèle standard.
- Le modèle cube `6+1` sur S³ à 60° donne un barycentre de `4/7`, mais ses distances possèdent **3 familles** au lieu des 6 familles de Polybrutus.

## Limite de physique cruciale : la jauge

Les différentes orientations du doublet Higgs sur la sphère de vide sont reliées par la symétrie de jauge. Sans espace-temps, connexion électrofaible et observable **invariante de jauge**, les 21 distances entre sept vecteurs de composantes du champ ne représentent pas des observables physiques indépendantes. Sept points sont ici **sept positions géométriques choisies**, pas sept nouveaux bosons.

## Dossier et provenance

- `RAPPORT_SCIENTIFIQUE.md` : équations, méthode, résultats exacts, contre-tests, limitations, prochaine expérience.
- `polybrutus_higgs_lab.py` : code complet sans dépendances, témoins intégrés et vérification PSD par Schur.
- `test_polybrutus_higgs_lab.py` : 10 tests unitaires, sans ajustement à une masse physique.
- `RESULTATS_2026_11_05.txt` et `TESTS_2026_11_05.log` : traces d'exécution.
- `PROVENANCE.md` : origine des matrices et des hypothèses.

**Ce laboratoire ne démontre pas une nouvelle particule, une nouvelle loi physique ni l'originalité dans la littérature.** Il montre ce qui est mesurable, ce qui est impossible avec certaines familles de fonctionnelles et quelle hypothèse supplémentaire il faudrait spécifier pour poursuivre.
