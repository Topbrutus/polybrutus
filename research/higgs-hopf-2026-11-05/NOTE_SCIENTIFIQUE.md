# Polybrutus et Higgs : projection de Hopf et limite de jauge (2026-11-05)

**Projet : Gabriel St-Pierre (Topbrutus). Recherche exacte non relue par les pairs. Aucune prédiction de boson ou de masse.**

## Construction
Les sept points unitaires de [SPECIFICATION.md](../../SPECIFICATION.md) sont des directions n_i sur S2. On peut les relever en sept doublets complexes normalisés psi_i dans S3, par la projection de Hopf n_i = psi_i^dagger sigma psi_i. La fibration S3 -> S2 est un résultat classique.

**21 paires :** |<psi_i|psi_j>|^2 = (1+n_i.n_j)/2 = 1-d_ij^2/4. Le multiensemble des probabilités est 3/4 (6 fois), 9/16 (2), 7/16 (3), 1/4 (6), 3/16 (2), 0 (2). Les paires antipodales sont (P1,P5) et (P2,P7).

**35 triangles :** B_ijk = <psi_i|psi_j><psi_j|psi_k><psi_k|psi_i> = (1+n_i.n_j+n_j.n_k+n_k.n_i + i det[n_i,n_j,n_k])/4. Il y a 25 valeurs non nulles et 10 nulles. Pour les valeurs nulles, la phase n'est pas definie. B_123 = 15/32 + 3 sqrt(3) i/32.

**Moment exact :** M = sum_i n_i n_i^T = diag(21/8,9/8,13/4), det M = 2457/256. Cauchy-Binet implique sum_{i<j<k} (Im B_ijk)^2 = det M / 16 = **2457/4096**. L'identite generale est connue; sa valeur pour Polybrutus est calculee ici.

## Contre-test d'invariance de jauge
Les B_ijk sont invariants sous rephasages U(1) independants et SU(2) globale; ils **ne** sont **pas** invariants sous SU(2) locale independante. Pour chaque psi_i=(a_i,b_i), G_i=[[a_i,-conj(b_i)],[b_i,conj(a_i)]] est dans SU(2) et G_i e1=psi_i. Les transformations locales G_i^dagger alignent toutes les orientations sur e1 et rendent chaque produit brut B=1.

Avec des liens U_ij se transformant comme h_i U_ij h_j^dagger, les produits B^U_ijk=(psi_i^dagger U_ij psi_j)(psi_j^dagger U_jk psi_k)(psi_k^dagger U_ki psi_i) sont invariants sous la transformation conjointe. Avec les champs fixes, les connexions plates U_ij=I et U_ij=G_i G_j^dagger donnent des valeurs distinctes. Elles sont equivalentes comme connexions seules, mais pas comme couples (champs fixes, connexion).

## Resultat et limite
Un laboratoire SymPy 1.14.0 local contient 19 tests PASS, un certificat minimal autonome, des resultats JSON et un manifeste SHA-256 verifie. Ces fichiers sont dans une archive fournie dans la conversation; **ce commit ne publie que cette note**.

La geométrie seule ne determine pas une observable electfaible : il manque une connexion SU(2)xU(1), des sites d'espace-temps et des couplages fixes independamment. Aucune prediction de la masse du Higgs ni d'un nouveau boson.

Prochaine action : fixer un modele de jauge electfaible independant et comparer des observables transportees a des temoins geometriques.

References : https://bohr.physics.berkeley.edu/classes/250/f15/hws/hw1.pdf ; https://www.sciencedirect.com/science/article/pii/S037596012600561X ; https://www-archive.ph.ed.ac.uk/sussp61/lectures/01_Sachradjda_StandardModel/standrews1.pdf
