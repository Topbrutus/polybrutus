# Note scientifique — boucle discrète à sept sites / doublet de Higgs

**Date :** 2026-11-05  
**Auteur de Polybrutus :** Gabriel St-Pierre (Topbrutus)  
**Type :** note de recherche expérimentale computationnelle ; revue par les pairs **non réalisée**.

## 1. Séparation indispensable des deux espaces

- **Espace de positions :** sept sites sur un anneau spatial `x_i`, reliés par `7` arêtes consécutives (pas les `21` paires). La distance spatiale des voisins est définie comme `a=1` dans l'exemple.
- **Espace interne :** à chaque site, une valeur `phi_i∈C²`, soit quatre composantes réelles, prise parmi les sept points d'une des configurations géométriques sur `S³`.

Aucun principe connu du modèle standard n'impose cette correspondance entre un point de l'anneau et un sommet Polybrutus. C'est la partie arbitraire la plus importante à tester.

## 2. Action et jauge

Pour les arêtes orientées `i→j` de l'anneau :

`E = (κ/a) Σ ||φ_i − U_ij φ_j||² + a Σ λ(φ_i†φ_i − v²/2)²`.

Transformation locale combinée :

`φ'_i=H_i φ_i`, `U'_ij=H_i U_ij H_j†`, avec `H_i=q_i R_i`, `|q_i|=1`, `R_i∈SU(2)`.

Alors `φ'_i−U'_ij φ'_j=H_i(φ_i−U_ij φ_j)`, donc l'action entière est **invariante de jauge** (la norme et le potentiel local sont conservés). Le produit ordonné des transports sur un cycle est transformé par conjugaison au point de base. Ainsi la trace de Wilson est elle aussi invariante.

On fixe `a=κ=λ=g=g'=1`, `v=√2`, `Y=1/2`; ces nombres sont purement illustratifs. Aucun champ de jauge continuum `A_mu,B_mu`, action Yang–Mills, ni équations temporelles n'est dérivé.

## 3. Calculs contrôlés

Dans le fond à liens identités, `E=Σ_{(i,j) voisins} d²_ij` pour les spinors normalisés de vide. Les valeurs exactes sont :

| Géométrie interne | Action de l'anneau, liens identités | Action, holonomie imposée `qR` |
|---|---:|---:|
| Polybrutus 3D | `31/2` | `743/50 − 12√3/25` |
| 4D A (cosphérique) | `18` | `418/25 − 32√2/75` |
| 4D C (cosphérique) | `16` | `384/25 + 108√183/1525` |
| Cube six faces + référence | `14` | `366/25 + 16√3/25` |

Le témoin 4D B est non cosphérique; on ne peut donc pas le placer directement dans l'hypersphère du vide unitaire en conservant toutes ses distances. Le programme conserve sa somme d'anneau **géométrique** `18` mais l'exclut des actions de Higgs à norme fixe.

La même holonomie a été utilisée pour **tous** les cas, sans optimiser les résultats. Elle est un paramètre externe du modèle-jouet, non une configuration physique mesurée.

## 4. Contre-tests qui empêchent la surinterprétation

1. **Permutation des sites :** en échangeant les orientations des sites 3 et 4, l'action identité de Polybrutus 3D change de `31/2` à `37/2`; les 21 distances non étiquetées demeurent les mêmes. La valeur dépend donc d'une identification supplémentaire entre points internes et sites spatiaux.
2. **Transformation de jauge locale appliquée seulement aux champs :** l'action change si les `U_ij` ne sont pas transformés aussi ; les distances brutes ne sont pas observables de jauge.
3. **Transformation simultanée champs + connexions :** les 7 énergies d'arêtes, l'action complète et la trace de Wilson sont égales avant/après sur chacun des quatre modèles sphériques (**contrôle exact**).
4. **Choix de la connexion :** des liens de transport `G_i G_j†`, construits *en fonction des champs* (et donc non admissibles comme prédiction indépendante), annulent toute énergie cinétique du vide. Cela montre à quel point la connexion influence le résultat.
5. **Fond à holonomie non triviale :** `tr(W)=18/25+24i/25`, alors que `tr(W)=2` pour les identités. Cette trace complexe est invariante pour la représentation choisie; elle ne constitue pas une courbure locale mesurée.
6. **Courbure radiale du potentiel :** `V(h)=h⁴/4+√2 h³+2h²`, `V''(0)=4` dans nos unités arbitraires. La relation générale standard `V''(0)=2λv²` est retrouvée; la valeur numérique est un choix, pas une prédiction du Higgs.

## 5. Portée et frontière physique

**Prouvé dans ce cadre :** existence d'une action discrète explicitement invariante sous transformations locales du doublet et des liaisons; calculs exacts pour quatre jeux de sept états; forte dépendance à la topologie, à l'appariement et au transport choisi.

**Non prouvé :** que Polybrutus soit imposé par le Higgs, qu'il existe sept états physiques distinguables, qu'un paramètre de masse du modèle standard soit déduit, que la topologie des sites soit celle de l'espace-temps réel, ou qu'une particule nouvelle soit prédite. Le modèle est un exemple `1D` statique, pas un calcul dans le modèle standard complet `SU(2)_L×U(1)_Y` couplé dynamiquement.

La symétrie de jauge peut rendre physiquement équivalentes plusieurs orientations internes ; une observable dépendant des arêtes exige donc un **champ de transport indépendant**. Sans justification indépendante des liens et du graphe, toute conclusion physique serait circulaire.

## 6. Condition d'abandon / prochaine expérience

**Prochaine action unique :** comparer ce même protocole à une cohorte de *témoins indépendants*, incluant des ensembles aléatoires de sept spinors unitaires, ainsi qu'à des permutations exhaustives de l'attribution des sept orientations aux sites, en conservant la connexion et les paramètres fixés. Établir une distribution de référence et une statistique annoncée à l'avance **avant** de chercher une interprétation physique.

Si aucune distinction robuste n'existe face aux témoins, **abandonner l'idée que les énergies de cet anneau révèlent une propriété physique distinctive de Polybrutus**. Si une distinction subsiste, demander une revue experte sur la pertinence du graphe et des transports avant de parler de physique.
