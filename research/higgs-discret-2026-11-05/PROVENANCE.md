# Provenance / Troisième œil — laboratoire discret sept sites

**2026-11-05 — lecture, intention, action, preuve, inconnues, reprise.**

## Sources relues

- `Topbrutus/polybrutus`, `main` lu au SHA `8e632d2c3c76da0dd85662edad155637617e7b00`.
- Branche `astra/polybrutus-higgs-variationnel-20261105` lue au SHA `2eaed69fa37c9cd8ef65c6ac8a36bd0d0f386fdd` avant cette étape ; PR #1 brouillon.
- Laboratoire variationnel: `research/higgs-variationnel-2026-11-05/polybrutus_higgs_lab.py`.
- Laboratoire de jauge: `research/higgs-jauge-2026-11-05/higgs_gauge_links.py`.
- `Topbrutus/topbrutus_one/memoria/gpt/sessions/ASTRA_POLYBRUTUS_HIGGS_JAUGE_2026-11-05.md`.
- `Topbrutus/topbrutus_one/memoria/THIRD_EYE_PROTOCOL.md`.

## Intention de Topbrutus

« continue » : passer du test algébrique de jauge à un modèle discret situé en espace-temps avec action covariante, paramètres fixés, tests et contre-exemples, sans toucher aux versions publiées et sans invention de résultat physique.

## Action

Code autonome (au sein du dépôt, avec deux modules sources frères), anneau spatial de sept sites, action avec connexions de `U(2)` dans la représentation doublet `SU(2)×U(1)_Y`, contrôle de jauge exact, Wilson non trivial, tableau de résultats, contrôle de permutation des sites, potentiel radial local. Statut du modèle : *démonstrateur mathématique*, **pas** la théorie électrofaible complète.

## Vérification et reproductibilité

`python -B higgs_lattice_7.py` puis `python -B -m unittest -v test_higgs_lattice_7.py`. Les résultats de la passe sont consignés dans les JSON et logs du même dossier, puis les empreintes dans `MANIFEST_SHA256.txt`. Voir les résultats et le ZIP indépendant pour preuve. Les SHA GitHub du laboratoire seront confirmés après publication sur la branche.

## Incertitudes

Aucun fond de jauge, aucun paramètre `a,κ,λ,v,g,g'`, aucune correspondance sept points ↔ sept sites n'est fixé par une mesure externe. Aucun test statistique contre une cohorte indépendante n'a encore été exécuté. La dimension physique de l'action n'est pas calibrée. L'anneau 1D ne décrit pas l'espace-temps physique 3+1D.

## Suite

Produire des témoins randomisés avec graine fixée, puis comparer toutes les `7!` permutations d'affectation sans changer les règles d'évaluation. Tant que ce contrôle n'est pas fait et que les observables ne sont pas motivées physiquement, **ne rien annoncer sur une nouvelle particule**.
