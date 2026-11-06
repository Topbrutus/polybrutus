# Polybrutus × Higgs — laboratoire des invariants de jauge

**Date :** 2026-11-05  
**Auteur et concepteur de Polybrutus :** Gabriel St-Pierre (Topbrutus)  
**Statut :** expérience mathématique exacte, **aucun boson ni mécanisme physique nouveau démontré**.

## Objectif

Tester une limitation fondamentale de l'analogie « 7 orientations de Polybrutus sur la sphère du vide de Higgs » : les orientations de norme constante d'un doublet complexe de Higgs sont reliées par les transformations locales `SU(2)_L` (et, plus généralement, par la jauge électrofaible `SU(2)_L × U(1)_Y`). Une distance ordinaire entre deux doublets situés en des sites distincts ne reste pas inchangée si l'on effectue indépendamment une transformation de jauge à chaque site.

Ce laboratoire introduit explicitement des **liaisons de jauge** `U_ij` et un invariant discret :

```text
phi'_i = h_i phi_i
U'_ij = h_i U_ij h_j†
I_ij = ||phi_i − U_ij phi_j||²
I'_ij = I_ij
```

**Attention :** les sept sites `i=0,...,6` sont une hypothèse du modèle jouet; ni leurs positions dans l'espace-temps ni les champs de jauge physiques ne sont déduits de Polybrutus.

## Reproduction

À partir du dépôt GitHub `Topbrutus/polybrutus`, branche de recherche, avec Python et SymPy 1.14.0 :

```bash
python -m pip install sympy==1.14.0
cd research/higgs-jauge-2026-11-05
python -B higgs_gauge_links.py
python -B -m unittest -v test_higgs_gauge_links.py
sha256sum -c MANIFEST_SHA256.txt
```

**Dépendance locale vérifiée :** `../higgs-variationnel-2026-11-05/polybrutus_higgs_lab.py`, figée dans la branche avant ce test. Aucun fichier du laboratoire antérieur n'est modifié.

## Étendue des contrôles

Quatre configurations : `Polybrutus_3D`, `4D_A`, `4D_C`, `Cube+1`. La configuration `4D_B` n'est pas cosphérique et n'appartient pas à l'hypersphère du vide à norme constante; elle reste un contrôle d'exclusion. Pour `4D_A` et `4D_C`, un plongement de `S³` est reconstruit **exactement** depuis les matrices de distances au carré (Cholesky, sphère circonscrite et vérification des 21 distances).

Le programme prouve la covariance du transport pour les 21 liaisons, ainsi que la platitude des 35 triangles de Wilson, par calcul matriciel explicite pour Polybrutus 3D et Cube+1, puis par lemme unitaire exact pour 4D_A/4D_C, avec un contrôle matriciel direct supplémentaire sur une arête de chacun.

Résultats complets : `RESULTATS_2026_11_05.json`. Tests : `TESTS_2026_11_05.log` (**17 réussis**). Discussion et limites : `NOTE_JAUGE_ET_HIGGS_2026_11_05.md`.

## Interprétation permise

- Les 21 distances originales restent mathématiquement exactes dans leurs géométries respectives.
- Il existe des transformations de jauge locales qui envoient les sept orientations normalisées vers un même vecteur `(1,0)`. Les distances *brutes* peuvent alors s'annuler.
- Les grandeurs construites avec des liaisons transformées simultanément restent invariantes : elles dépendent donc de **l'ensemble (champs, connexion)**, pas des seules orientations.
- L'énergie à deux corps avec liaisons identités conserve `195/4` pour les trois Polybrutus sphériques, et la somme triple transportée vaut respectivement `24435/64`, `23835/64`, `24273/64` pour 3D, A et C.
- Une autre connexion *plate* (transport adaptatif `U_ij=G_i G_j†`) annule toutes les différences covariantes, sans changer les sept orientations.

**Aucune de ces constatations ne prédit la masse du Higgs ni ne démontre qu'une structure à sept sites est imposée par le modèle standard.**

## Prochaine action unique

Définir **avant tout ajustement** : une topologie ou des positions d'espace-temps, une discrétisation de `|D_μΦ|²` avec ses couplages physiques et une configuration de jauge indépendamment fixée; puis comparer les mêmes observables invariantes pour Polybrutus et des témoins.
