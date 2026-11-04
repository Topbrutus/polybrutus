# POLYBRUTUS — Définition mathématique de référence

**Concepteur / auteur du modèle documenté : Gabriel St-Pierre** (Topbrutus). Version de travail 0.1.0. Géométrie **exacte**; antériorité mathématique **non établie**.

## 1. Espace et paramètres

Dans l'espace euclidien R³, choisir un rayon R strictement positif, un centre de référence O=(0,0,0) et un axe nord–sud Z. Les sept points appartiennent à la sphère S_R: x²+y²+z²=R².

## 2. Configuration analytique

Les coordonnées ci-dessous sont exactes. La deuxième, troisième et quatrième position ont une colatitude de 60° depuis +Z; la sixième et la septième, une colatitude de 60° depuis −Z. Les trois points supérieurs ont les azimuts 0°, 120°, 240°.

| Point | x | y | z |
|---|---|---|---|
| P1 Nord | 0 | 0 | R |
| P2 Haut-1 | √3 R/2 | 0 | R/2 |
| P3 Haut-2 | −√3 R/4 | 3R/4 | R/2 |
| P4 Haut-3 | −√3 R/4 | −3R/4 | R/2 |
| P5 Sud | 0 | 0 | −R |
| P6 Bas-1 | √3 R/2 | 0 | −R/2 |
| P7 Bas-2 | −√3 R/2 | 0 | −R/2 |

Le centre **O** est distinct du barycentre **G**. Ni O ni G n'est un huitième sommet.

## 3. Graphe de liaisons

Chaque paire non ordonnée de sommets définit une liaison. Il y a C(7,2)=21 segments. Tous les segments sont étudiés ; ils ne sont **pas** tous des arêtes de la frontière d'un polyèdre fermé. Cette construction ne doit donc pas être présentée sans autre preuve comme un solide plat à sept faces.

## 4. Identité des rayons

En écrivant x=a√3 R, y=bR et z=cR avec a,b,c rationnels, on vérifie pour chaque sommet 3a²+b²+c²=1. Cela prouve exactement que les sept sommets appartiennent à la même sphère.

## 5. Invariants des six distances

Pour chaque paire i<j, le carré de distance est d²=(x_i−x_j)²+(y_i−y_j)²+(z_i−z_j)². Poser k=4d²/R².

L'ensemble exact est **{4,7,9,12,13,16}**, avec les multiplicités **{4:6, 7:2, 9:3, 12:6, 13:2, 16:2}**. La somme des multiplicités donne **21**.

Les six longueurs sont, dans l'ordre : R ; R√7/2 ; 3R/2 ; R√3 ; R√13/2 ; 2R.

## 6. Barycentre (non confondu avec le centre de la sphère)

ΣP_i=(0,0,R/2). Ainsi **G=(1/7)ΣP_i=(0,0,R/14)**. La relation G_z/R=1/14 est un invariant de changement d'échelle.

## 7. Conditions de reproduction

Le dossier est validé mathématiquement si le calcul indépendant vérifie à chaque exécution :
1. sept sommets distincts ;
2. rayon commun R ;
3. 21 liaisons ;
4. famille des k exacte ;
5. multiplicités exactes ;
6. barycentre exact ;
7. colatitude 60° pour les cinq sommets intermédiaires et alignement des deux points du bas dans XZ.

Exécuter `python verify_polybrutus.py`. Les composantes de x utilisent le symbole √3 mais les normes carrées sont calculées exclusivement avec des nombres rationnels exacts, sans arrondis.

## 8. Sur le mot « parfait » et sur la nouveauté

Le mode « exact 60° » impose une **construction analytique** et des **relations internes vérifiables**. Il n'existe dans ce dépôt **aucune preuve** qu'elle soit optimale pour un objectif universel non défini, ni qu'elle soit inédite dans la littérature. Le nom propre **Polybrutus** identifie cette réalisation documentée et n'établit pas, par lui-même, un droit exclusif sur un objet mathématique.

## 9. Référence et provenance

Document construit à partir des fichiers de travail **TOPBRUTUS_7POINT_SPEC.md** et **generateur_7_pics_snap_exact_60_invariants.py**, retrouvés dans la bibliothèque des travaux de Gabriel St-Pierre. Voir [PROVENANCE.md](PROVENANCE.md).
