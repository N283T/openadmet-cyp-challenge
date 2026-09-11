# Data and prediction integrity audit

The checks below found no molecule/label permutation, stale regression data,
missing inverse target transform, or submission packaging error. They do not
establish the cause of the leaderboard gap.

## Official data freshness

Run `20260911T141511080701Z-fd93c411` resolved the live Hugging Face dataset to
`3ac9c5dbb83eec5780ec7fa511908698cfe1396d` and downloaded all five files at that
revision, preserving the original local raw files. Test, inhibition, TDI and
single-concentration files are byte-identical to the local files.

Emax changed from 6,146 to 6,145 rows. Run
`20260911T141633727753Z-44127598` compared full rows with null-aware matching:
the only removal is `OCNT-2329004`, SMILES `N.N.[Cl][Pt][Cl]` (cisplatin), with
no additions or edited surviving rows. This exact structure was already
quarantined by our ChEMBL standardization. It did not enter model features.
Consequently this update does not change the current regression training data.

Source: https://huggingface.co/datasets/openadmet/cyp-challenge-train-test/tree/3ac9c5dbb83eec5780ec7fa511908698cfe1396d

## Molecular identity

Recomputed standardization for all 6,897 unique raw SMILES. Stored canonical
structures, compound hashes and statuses match exactly: 6,882 unchanged, 14
changed, one excluded. All 6,896 retained compounds match the saved compound
table. No molecule name maps to multiple standardized compound IDs.

The 14 changes comprise ten carboxylate neutralizations and four sulfoxide
representation normalizations. None loses heavy atoms; detected assigned
stereocenter counts are retained. Raw/standardized strings and counts are saved
for inspection. This is a consistency and structural-change audit, not a claim
that every raw assay structure has been independently authenticated chemically.

## Labels and submitted predictions

Reconstructed all 6,525 direct-inhibition pIC50 labels and confidence bounds from
raw training CSVs through molecule-name/SMILES mappings. They exactly match
the measurement table. All 6,384 current TabICL OOF labels and bounds also
match those raw values. No additional IC50-to-pIC50 conversion is applied.

All 750 submitted molecule names and raw SMILES match the official test file.
Recomputed the keyed five-fold means from 15,000 saved test prediction records;
they match the assembled predictions within 1e-12. Packaging changes no values.
The actual prepared submission and packaging output have the same SHA-256:
`6a37411ff8a7c6a28a29f76752c2b9bc790e9fbc9766a836ba3ec5354894576d`.

## Feature and model replay

Frozen embeddings have complete, unique compound IDs. Recomputed 142 structures
(128 seeded random compounds plus all 14 changed structures) from SMILES and
the saved CheMeleon checkpoint. Maximum difference from stored embeddings:
4.17e-7.

Reloaded the adopted pair-0/1 auxiliary best checkpoint and recomputed the same
142 compounds. Head order is CYP1A2, CYP2C9, CYP2D6, CYP3A4. Applied the saved
training standard deviations and means exactly once. Maximum difference from
exported log2fc predictions: 9.54e-7. Both checks pass atol/rtol 2e-5.

Source inspection confirms downstream fitting uses raw pIC50 labels, shared
training-fitted PCA/imputation for validation and test, and identical concatenation
order. The existing replacement audit already checks reconstructed features
against saved LGBM predictions and auxiliary excluded-fold memberships.

This new replay samples one of ten auxiliary models; it does not numerically
replay every auxiliary checkpoint, every Mordred descriptor, or all downstream
TabICL fitted models. No new training or submission was performed. Remaining
modeling differences, including the full feature recipe and reduced training
sets, are not explained away by these passing integrity checks.
