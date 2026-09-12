# Ligand-only geometry and charge proxy experiment

The report suggested CYP1A2 aromatic planarity and CYP2D6 basic-nitrogen/aromatic
spatial relationships. This experiment isolates those signals from the existing
frozen Jeremy 600-dimensional embedding plus four nested primary-screen log2fc
predictions. It is separate from Jazzy and uses neither protein structures nor docking.

## Feature construction

Feature run: `20260912T092823285236Z-3bc4bd9d`.
All 6,896 standardized compounds are retained. RDKit 2026.03.6 ETKDGv3 generates
three explicit-hydrogen conformers (seeds 11, 42, 89), with no force-field
minimization. 20,675 of 20,688 conformer attempts succeeded; 13 failed embedding.
The run completed in approximately 34 seconds. Seed-level descriptors and failures
are saved in `geometry-seeds.parquet`; aggregated descriptors are keyed by compound ID
in `geometry-features.parquet`.

Ten topology/charge descriptors comprise formal charge, positively charged nitrogen
count, basic-N structural proxy count, aromatic N/atom/ring counts, rotatable bonds,
minimum/maximum nitrogen Gasteiger charge and total absolute heavy-atom Gasteiger charge.
The structural basic-N proxy selects aliphatic nitrogen excluding N attached to
C/S/P double-bonded to O/S/N, plus pyridine-like aromatic nitrogen. This deliberately
simple pattern is incomplete (for example, it can omit amidine/guanidine nitrogens);
it does not estimate pKa or enumerate pH-dependent microstates.
Gasteiger partial charge likewise is not a pKa prediction.

Thirteen conformer descriptors comprise whole-molecule plane RMS/max deviation,
NPR1/NPR2, radius of gyration, asphericity, eccentricity, global aromatic plane RMS
deviation, maximum aromatic ring centroid distance, mean absolute ring-normal dot
product, minimum/mean basic-N to aromatic-centroid distance (excluding an N's own
ring), and basic-N inverse nonbonded heavy-atom crowding within 4 angstroms. The
last descriptor is an exposure proxy, not solvent-accessible surface area.
Both means and population standard deviations across successful seeds are saved.
Inapplicable descriptors, such as nitrogen charge for nitrogen-free molecules,
remain missing, as do failed conformers; compounds are never removed.

## Controlled evaluation

Readout run: `20260912T092909986920Z-c5d42521`.
Four prespecified arms use identical LGBM L1 settings (300 trees, 15 leaves,
minimum leaf size 20, learning rate 0.05, L2=1, seed 42+fold):
control; ten topology/charge proxies; those proxies plus 13 geometry means;
and those proxies plus means and 13 conformer standard deviations.

All arms use the fixed parent assignments `a3923255`, frozen embeddings `c060d197`,
and updated pair-excluded primary predictions `d272774c`. Nonconstant feature
selection and median imputation use outer training rows only, with missing indicators
appended for retained columns. Anchors remain training rows. The control must
reproduce `1e7001de` bit for bit. Each arm emits 6,384 keyed OOF predictions with
endpoint, fold, split, labels and confidence bounds; fitted models, preprocessing
states and training memberships are retained.

## LGBM results

Independent audit: `20260912T093150407906Z-615e167c`. Control reproduction was
bitwise exact; all 6,384 keys, labels, intervals and validation memberships passed.

| Arm | 1A2 ST-RAE | 2C9 | 2D6 | 3A4 | Macro |
|---|---:|---:|---:|---:|---:|
| Control | 0.691998 | 0.494360 | 0.842287 | 0.442934 | 0.617895 |
| Topology/charge proxies | 0.691492 | 0.489857 | 0.847586 | 0.441300 | 0.617559 |
| Proxies + geometry means | 0.685658 | 0.484669 | 0.844729 | 0.438636 | 0.613423 |
| Proxies + geometry means/std | 0.689052 | 0.482692 | 0.839852 | 0.438107 | 0.612426 |

These inexpensive descriptors have a small standalone effect. Means help 1A2;
conformer variability helps 2D6 relative to the same control, while topology-only
charge proxies do not resolve 2D6. This is a modest feature-family result rather
than evidence that a particular descriptor captures the biological mechanism.

The fixed pre-Caruana `ae562a59` ensemble has macro ST-RAE 0.586931. Adding each
geometry arm at prespecified 25% or 50% does not improve its four-endpoint macro:
best is 0.588380 for 25% mean/std. Focus endpoint gains are tiny: 25% geometry means
changes 1A2 from 0.670809 to 0.669463; 25% mean/std changes 2D6 from 0.822681 to
0.822186. These separately observed endpoint choices are not a newly validated
combined model. The reference here is explicitly the earlier fixed ensemble, not
the more recent global Caruana crossfit result. All candidate OOF vectors are
available for the common global comparison after the report experiment round.

Feature quality audit `20260912T093535838085Z-f2107add` confirmed there are no
remaining disconnected standardized structures. 6,890 molecules have all three
successful conformers, two have two, one has one, and three have none. The last
three remain in every applicable model through training-only imputation and
missing indicators. Per-column missing/nonfinite counts are saved in `quality.json`.

## TabPFN confirmation

Paired run `20260912T093227551001Z-ed3e4d03`, independent audit
`20260912T093909318053Z-476204b3`: the same control versus all 36 descriptors,
TabPFN 2.6 with eight estimators, same fold seeds and preprocessing. The original
PFN control `568bc526` also reproduced bit for bit (maximum difference 0).
No hyperparameter search was performed.

| Arm | 1A2 ST-RAE | 2C9 | 2D6 | 3A4 | Macro |
|---|---:|---:|---:|---:|---:|
| TabPFN control | 0.705017 | 0.495867 | 0.874374 | 0.439395 | 0.628663 |
| TabPFN + proxies/geometry means/std | 0.692620 | 0.485648 | 0.877186 | 0.436212 | 0.622916 |

The PFN result supports a small 1A2 signal, but does not improve 2D6. A fixed 25%
addition to `ae562a59` changes 1A2 to 0.669797, while the overall macro worsens to
0.590082. The geometry PFN arm improves neighbor-subset macro when added at 25%
(0.675304 to 0.674066), but that alone is insufficient for adoption.

Retain the three LGBM geometry variants and the one PFN geometry variant as
additional global ensemble candidates. There is no large focus-endpoint gain,
no newly submitted model, and no evidence here that ligand geometry fully captures
protein binding, protonation or assay-specific effects. Check/formats for the six
new `report_geometry*.py` scripts pass Ruff.
