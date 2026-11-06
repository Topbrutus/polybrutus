# Provenance et intégrité — 2026-11-05

- Définition initiale : `Topbrutus/polybrutus` / `SPECIFICATION.md` (branche `main`), 7 points de R³ sur S².
- Source directe des doublets de S³ : `research/higgs-jauge-2026-11-05/higgs_gauge_links.py` (branche `astra/polybrutus-higgs-variationnel-20261105`). SHA-1 du **blob GitHub** vérifié `981ae8f7c74123f94f0590c5b41e581c7d5ec291` ; SHA-256 local `c33952388df8e45bce9afd82b3fe7aa125d115b3f33e5608eb959df46b4f9e1d`.
- `GEOMETRIES_EXACTES.json` recopie **les expressions symboliques exactes**, sans les recalculer ni inventer de quatrième dimension; les contrôles SymPy valident normes et 21 distances.
- Le modèle de jauge précédent existant `research/higgs-discret-2026-11-05/` avait déjà une action statique SU(2)×U(1), des fonds fixes et des témoins. Le présent dossier **ne revendique pas la première construction**; il ajoute l'analyse spectrale du vide neutre, un minimum et une borne inférieure exacts, ainsi que des tests d'identifiabilité sur les quatre témoins.
- Paramètres du nouveau protocole : `q=(3+4i)/5`, 4 fonds, 7 sites, action dimensionnelle fictive avec `a=kappa=lambda=1`. Tous sont des choix mathématiques et non des valeurs issues de données physiques.
- Le titre « Polybrutus × Higgs » indique l'objet étudié, et non une affirmation de lien expérimental.
- Référence sur la théorie établie : Weinberg–Salam / mécanisme BEH et groupe électrofaible SU(2)L×U(1)Y ; formalisme standard de réseaux de jauge. Le résultat spectrale propre à cet anneau se vérifie directement avec les formules et le code fournis, sans recours à une publication non vérifiée.
- Tests : suite `unittest`, logs, JSON, archive vérifiée dans un dossier vierge ; contrôler le manifeste SHA-256 avant toute réutilisation.
- Aucune écriture requise dans `main` de Polybrutus ou dans Zenodo. Ce dossier peut être ajouté **uniquement** à la branche de PR en brouillon si la version et les tests sont vérifiés.
