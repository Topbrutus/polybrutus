# Polybrutus × Higgs — modèle électrofaible statique à sept sites (2026-11-05)

**Responsable du projet :** Gabriel St-Pierre (Topbrutus). **Statut :** laboratoire mathématique exploratoire, ni théorie électrofaible complète ni résultat évalué par les pairs. Les calculs numériques sont dans des unités arbitraires.

## Le problème précis

Les orientations des sept points de Polybrutus peuvent être codées dans sept doublets complexes normalisés. Cette analogie détermine-t-elle une observable électrofaible **sans** ajuster la connexion de jauge après avoir examiné les points ? Le contrôle négatif consiste à autoriser, séparément et explicitement, une connexion plate adaptée aux champs : cela doit effacer toute différence.

## Protocole fixé avant les résultats

- Sept **sites physiques fictifs** sur un anneau unidimensionnel périodique, ordonnés 0→1→…→6→0, pas de maille `a=1`.
- Aux sites, un doublet complexe `phi_i=(phi_i^+,phi_i^0)`, `||phi_i||=1` (vide imposé), `lambda=1`, `kappa=1`; potentiel `sum_i (||phi_i||²−1)²`.
- Sur chaque arête orientée, une matrice SU(2) `S_ij` et une phase U(1)_Y agissant sur le doublet `z_ij`.
- Action statique `E = sum_edges ||phi_i−z_ij S_ij phi_j||² + sum_i (||phi_i||²−1)²`.
- Pour changement de jauge local indépendant aux sites : `phi_i'=t_i h_i phi_i`, `S_ij'=h_i S_ij h_j†`, `z_ij'=t_i z_ij conj(t_j)`. Ici la phase `t_i` désigne directement la représentation U(1)_Y sur le doublet, avec hypercharge conventionnelle `Y=+1/2`.
- Quatre connexions prédéfinies, toujours identité sur les six premières arêtes. Sur l'arête 6→0 :

| Fond | Facteur SU(2) | Facteur U(1) sur le doublet | Critère |
|---|---|---|---|
| `identity` | `I` | `1` | référence |
| `su2_only` | `diag(q, q*)` | `1` | holonomie générique qui n'a aucune valeur propre 1 |
| `neutral_em` | `diag(q, q*)` | `q` | holonomie combinée `diag(q²,1)` laissant la composante neutre inchangée |
| `hyper_only` | `I` | `q` | contrôle hypercharge seulement |

`q=(3+4i)/5` est fixé pour permettre des identités exactes avec rationnels. **Ce n'est pas une valeur expérimentale de couplage électrofaible.**

- Jeux géométriques pris dans le laboratoire **antérieur** : Polybrutus 3D plongé dans S³; témoins cosphériques 4D A et 4D C avec exactement les mêmes 21 distances non étiquetées; témoin cube + point. Le témoin 4D B **non cosphérique** est écarté de ce protocole de norme fixe, pas de la recherche géométrique.
- Paramètres et connexions fixés indépendamment des quatre jeux. La correspondance points↔sites n'a pas de justification physique et reste un contrôle obligatoire.

## Résultat discriminant et limite

- Le modèle donne une holonomie SU(2) non triviale avec un Higgs neutre covariamment constant si la phase U(1)_Y compense la phase SU(2), ce qui reproduit la structure de l'U(1) électromagnétique non brisé.
- Holonomie SU(2) seule : minimum **strictement positif** parmi toutes les sept orientations unitaires. Il est calculé par le spectre exact du laplacien tordu.
- Si les liens sont librement ajustés aux sept champs, **tous** les jeux admettent une connexion SU(2) plate de coût nul, donc le modèle ne sélectionne pas Polybrutus.
- Avec un fond fixe, les énergies dépendent de l'attribution des sept doublets aux sept sites, qui n'a pas été déduite de la physique.

**Aucune identification avec un boson réel, aucune masse du Higgs prédite, aucune équation de mouvement 3+1D testée.** Le résultat principal est une contrainte de cohérence et d'identifiabilité d'un modèle-jouet, pas une découverte de physique.

## Reproduction (autonome)

Python 3.12+ et `numpy`, `sympy` ; versions utilisées dans la session : NumPy et SymPy disponibles localement (SymPy 1.14.0). Aucun téléchargement ni fichier d'un laboratoire voisin n'est requis à l'exécution. Les géométries exactes sont recopiées et sourcées dans `GEOMETRIES_EXACTES.json`.

```bash
python -m pip install numpy sympy
python -m unittest -v test_electrofaible_sept_sites
python electrofaible_sept_sites.py > verif_resultats.json
```

Vérifier `RESULTATS_2026_11_05.json`, `TESTS_2026_11_05.log`, `MANIFEST_SHA256.txt` et `NOTE_SCIENTIFIQUE.md`. Tous les 5040 ordres des sept champs sont parcourus dans l'expérience de permutations identité. La fonction `run()` réévalue les invariances sous transformations locales pour les 4 configurations × 4 fonds.

**Références internes :** [Polybrutus](https://github.com/Topbrutus/polybrutus/blob/main/SPECIFICATION.md) ; [jauge et transports](https://github.com/Topbrutus/polybrutus/tree/astra/polybrutus-higgs-variationnel-20261105/research/higgs-jauge-2026-11-05) ; [laboratoire statique précédent](https://github.com/Topbrutus/polybrutus/tree/astra/polybrutus-higgs-variationnel-20261105/research/higgs-discret-2026-11-05).
