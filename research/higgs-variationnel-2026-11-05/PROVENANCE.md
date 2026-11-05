# Provenance et portée des preuves — 2026-11-05

## Données externes contrôlées

- `Topbrutus/polybrutus`, branche `main`, HEAD lu : `8e632d2c3c76da0dd85662edad155637617e7b00`.
- `SPECIFICATION.md` : données originales P1..P7 et six familles de distances, version inspectée SHA blob `765c2d15c4a7f063ce506da033bd4f3a8049ef3b`.
- `research/distance-reconstruction-3d-4d/NOTE_RELECTURE_SCIENTIFIQUE_2026_11_05.md` : témoin A (coordonnées Q1..Q7) et statut non évalué par les pairs.
- Artefact local antérieur : `/mnt/data/polybrutus_certificat_trois_classes_4d_2026_11_05.py` : matrices B/C et preuve exactes de Gram, sphère et rang. Les matrices sont intégrées au vérificateur actuel et recalculées indépendamment; le balayage aléatoire de recherche n'est **pas** relancé ici.
- Artefact local antérieur : `NOTE_CUBE_SEPT_S3_HIGGS_2026_11_05.md` et `verification_cube_sept_higgs_s3_2026_11_05.py` : sept points cubiques sur S³, angle imposé 60°.

## Ce qui est réellement vérifié par ce dossier

Le programme reconstruit les matrices de distances à partir des coordonnées Poly3 et A, charge les matrices B/C en entiers, contrôle les 21 distances avec multiplicités, la PSD et le rang via complément de Schur exact, les sphères, les énergies pair/triple et les gradients projetés. Les 10 tests locaux ont été exécutés, et le log d'exécution est inclus.

## Ce qui est exclu

- Aucune expérience de physique des particules.
- Aucune équivalence démontrée entre les sept points et le champ de Higgs.
- Aucune preuve de priorité, de publication scientifique évaluée par les pairs, d'exhaustivité des réalisations R⁴.
- Aucune manipulation des DOI, de Zenodo ni des fichiers de référence du boson de Brutus.
- Aucune observation de serveur Antmux ou de simulateur distant; travail local en lecture seule sur les sources.

## Prochaine action unique

Établir un observable qui résiste à une transformation de jauge indépendante aux points du champ, ou conclure explicitement qu'aucune interprétation physique n'a été produite. Ne pas paramétrer un modèle à partir de la masse du Higgs observée.
