# Exact community CYP data reconstruction

The public author's data and historical18-head union can be reconstructed directly. The resulting union has **31,670 rows and161,400 labels**, exactly matching the published working notes. ChEMBL37 is supplied by the author's actual published table, so no ChEMBL36 fallback or approximate re-query is used.

## Immutable sources

- Source/history snapshot: `20260912T103808283993Z-48050646`, repository commit `a1d1b46e290e62014ee3f0dfde5667eb5e93282b`.
- Historical18-head model/producer: commit `971f6a7f24d86e59bab219347c3e14c006ae4e35` (18-head definition remains in the subsequent4cf07a8c formatting revision).
- Actual author public CSV download: `20260912T103958592504Z-a85b8d52`, including HTTP ETag, Last-Modified and full source URLs in the manifest. Run inventory preserves file hashes.
- Polars reconstruction and canonical-ID bridge: `20260912T104336148000Z-ae7e37fc`.
- Author blinded750 structure-only CSV: `20260912T104508219654Z-443b8901`.
- Independent blinded identity/SMILES comparison: `20260912T104624702551Z-5cc21663`.

The [public data client](https://github.com/SuperCowPowers/workbench/blob/a1d1b46e290e62014ee3f0dfde5667eb5e93282b/src/workbench/public_data.py) explicitly reads the anonymous `workbench-public-data` S3 bucket. The downloaded tables are the exact datasets requested by the author's producers, not scraped summary statistics.

## Four member contracts

Native target names and author-supplied SMILES are retained. Canonical project compound IDs are attached to challenge rows via exact molecule-name aliases. Public rows retain their author identity and get stable separate IDs. Their connectivity remains available for outer-fold exclusion even when a public row corresponds to challenge chemistry outside a member's original challenge cohort.

| Member table | Native head contract | Training rows with any selected label |
|---|---|---:|
| `challenge8.parquet` | Four direct pIC50 + four log2fc | 4,905 |
| `union18.parquet` | Same8 + five ChEMBL37 potency + five qHTS maximum-response | 31,670 |
| `d6_isoform6.parquet` | D6 direct pIC50,log2fc,TDI pIC50,**direct-condition Emax**,ChEMBL potency,qHTS maximum response | 23,264 |
| `d6_single1.parquet` | D6 direct pIC50 alone | 1,493 |

Head ordering is serialized in `contract.json`. Each table includes original CI columns for evaluation, but those are not model features. Exact model-specific source cohorts are also saved as `*-cohort.parquet`.

The18-head union does not contain TDI or Emax heads and does not contain Veith potency heads. Those appeared in later script versions. The D6 six-head member uses the later side-arm producer: TDI is the side-table base, Emax joins by name, and1,240 TDI-only molecules are appended to the inhibition cohort. The single Emax-only molecule is excluded exactly as the source code specifies. Its challenge cohort is6,145 and its ChEMBL/qHTS union cohort32,877; selecting any of the six D6 labels leaves23,264 training rows. The script's approximate16k prose is stale relative to its actual explicit `notna().any(axis=1)` selection and currently published data.

## Reproduced coverage

The public inhibition CSV has4,905 rows; the single-concentration CSV has17,504 rows across4,376 compounds. The source producer **left joins onto the inhibition cohort**, retaining4,375 measured compounds and17,500 log2fc labels. Our previous broader auxiliary preparation retained the extra single-concentration molecule. That difference is now removed for the historical8/18-head models.

| Target family | 3A4 | 2C9 | 2D6 | 1A2 | 2C19 |
|---|---:|---:|---:|---:|---:|
| Direct pIC50 | 2,335 | 1,285 | 1,493 | 1,412 | — |
| log2fc | 4,375 | 4,375 | 4,375 | 4,375 | — |
| ChEMBL potency | 15,777 | 10,554 | 8,415 | 10,968 | 10,991 |
| qHTS maximum response | 16,134 | 16,134 | 16,134 | 16,134 | 16,134 |

Total:6,525 direct +17,500 log2fc +137,375 public =161,400 labels. The D6 six-head model additionally has1,497 TDI pIC50 and1,493 direct-condition Emax labels.

## Faithful preprocessing, including public overlaps

The downloaded ChEMBL37 table contains25,841 rows. The source keeps the first record per supplied InChIKey connectivity block, yielding24,918 skeletons. It does not reaggregate values across those connectivity duplicates. The upstream loader uses exact nonnull pChEMBL values with null data-validity comment, potential_duplicate=0 and confidence>=8. Unlike our original stricter IC50-only ChEMBL36 export, it does not exclude imported qHTS records or restrict all values to direct-inhibition IC50. The author's column name `pic50_chembl` is preserved as published; this does not reclassify all underlying potency types as literal IC50.

Veith's85,535-row published table is averaged per exact author SMILES and isoform (matching pandas pivot_table's mean), then clipped to[-150,50]. SMILES rows are lexicographically ordered as in the source pivot. Connectivity is calculated on the author SMILES, and the first row per connectivity is retained. This yields16,131 unique skeletons. The union's16,134 labels per isoform include duplicated challenge connectivity rows, matching the original join behavior; they are not invented extra public molecules.

Source overlaps are merged as additional heads on existing challenge rows, not discarded:

| Public source | Inhibition-cohort overlaps | All training-compound overlaps | Blinded-compound overlaps |
|---|---:|---:|---:|
| ChEMBL37 | 185 | 217 | 0 |
| Veith | 160 | 183 | 1 |

ChEMBL's185/0 counts reproduce the blog's claim. The one Veith overlap carries a public historical assay response, not a hidden challenge label, and is retained to match the source. No hidden blinded pIC50 is downloaded or inferred. For CV, the trainer must remove **all head labels for every outer-held-out connectivity**, including public aliases and shared rows. `challenge-connectivity-folds.parquet` supplies the complete6,896-molecule bridge.

## Validation and prediction identities

Every one of the6,525 source direct labels and both CI bounds matches the current canonical project measurements within1e-9. Thus all1,457 non-anchor D6 OOF labels remain available; no validation molecule is dropped to force a source count.

`prediction-compounds.parquet` contains all6,896 canonical IDs. It uses author challenge SMILES where available and original official structures otherwise. A separate download of the author's actual blinded CSV confirms750 exact test IDs and **zero SMILES changes** relative to this prediction table. The blinded CSV has only molecule names and SMILES.

The trainer receives these prepared tables directly. No descriptor additions, transfer-learning initialization, target collapsing or alternate data weighting is introduced by this preparation. Training/runtime fidelity is documented separately by the reproduction trainer audit.

## Files

- `scripts/community_source_freeze.py`: source commit/history and public loader snapshot.
- `scripts/community_public_fetch.py`: exact published tables, response metadata, historical producer.
- `scripts/community_prepare.py`: source-faithful joins, four cohorts/contracts, canonical bridge and label checks.
- `scripts/community_prediction_audit.py`: independent blinded identity/SMILES verification.
- `configs/community_reproduction/`: immutable-run configuration sources.

All four scripts pass Ruff. Source data, generated tables and runtime artifacts remain outside Git.

## Independent native trainer data and mask audit

Recorded audit `shared/20260912T105536997597Z-08a9df69` reads the actual source tables, effective trainer configurations, emitted assignments, and available per-fold membership files as declared inputs. It confirms the selected head order and every finite-label count against the data contract, and reconstructs all five train/validation masks independently from saved assignments.

| Native run | Member | Source rows | Connectivity groups crossing folds |
|---|---|---:|---:|
| `fec30715` | Single direct, parent CV | 1,493 | 0 |
| `423d2caf` | Challenge auxiliary8, parent CV | 4,905 | 0 |
| `288cae8b` | Union18, parent CV | 31,670 | 0 |
| `e6e93fbf` | D6 isoform6, parent CV | 23,264 | 0 |
| `03b86209` | Single direct, source scaffold CV | 1,493 | 2 |

For all four parent-CV members, no training row shares connectivity with a held-out challenge molecule in any fold, including public aliases and additional public heads. Anchors remain in training. Every expected held-out challenge identity present in each member's input is in its assigned validation fold. The native trainer excludes complete rows, so this exclusion covers every observed auxiliary label, not only the scored direct head.

The source-scaffold single control retains the original splitter's assignments. Two connectivity groups cross its folds: one train/validation overlap in fold0, one in fold1, and two in fold4. These actual-source conflicts are exported as `20260912T105256414096Z-03b86209-connectivity-conflicts.parquet`. The source control should therefore be identified separately from the clean parent-CV comparison; its behavior was not silently modified to make the source audit pass. The first audit `bc769d01` failed its strict zero-cross-connectivity assertion on this source control and remains recorded; the second audit completes the diagnostics for all five runs.

At the audit snapshot, actual membership artifacts were available for all five folds of both single models, folds0–2 of auxiliary8, and fold0 of union18 and isoform6. Every available artifact agrees exactly with its saved assignment mask. All five masks for every model were checked from complete assignments; later memberships and final prediction coverage remain the trainer/output audit's responsibility.

The native trainer initially used shortened endpoint names when exporting scored OOF, producing empty scored exports despite valid preserved all-head predictions. The trainer owner identified this export-only mismatch and is recovering canonical endpoint OOF in a separate recorded run. This data audit does not interpret those initial empty exports as a model-training or label-coverage failure.

The four member training runs use source-native architecture, loss and fitting code with explicitly adapted parent CV. The separate single control uses source-native scaffold CV. Reproducing the data and training code does not make these two CV policies interchangeable. Metric-denominator comparisons, live leaderboard verification and eventual full-data refits are documented by the integration work rather than inferred from this data audit.
