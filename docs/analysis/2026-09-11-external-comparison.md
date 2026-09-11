# First matched external-data regression comparison

External auxiliary training improves the point estimate of pooled UMAP MA-ST-RAE
from **0.824695 to 0.797301**, but remains worse than the adopted TabPFN 2.6
baseline (**0.724704**). Keep the external dataset; do not replace the baseline
with this first direct-neural recipe.

The experiment tests one jointly trained auxiliary recipe. It does not establish
that the external dataset is ineffective for other training schedules or readouts.
The pooled improvement is not unambiguously resolved by fixed-OOF cluster
resampling, and only one training seed per fold was used.

## Method and evidence

The [external-data report](../external-data.md) specifies original source retrieval,
curation, decontamination, loss weights, stopping, and limitations. Both arms use
CheMeleon initialization and an identical 18-head architecture. The control fits
four challenge regression and four log2fc objectives; external heads are inactive.
The augmented arm additionally fits five ChEMBL IC50 and five PubChem response
objectives on challenge-disjoint chemistry. Classified TDI labels are not used.

Membership is identical between matched arms for every fold. Outer assignments
remain fixed UMAP; connectivity groups shared with an outer holdout are removed
from fitting, and inner validation uses disjoint connectivity groups. The external
dataset excludes all challenge connectivity, including test and auxiliary-only
molecules. Stopping depends only on inner-validation direct regression MSE.

Audit: `20260911T055219670967Z-ba1666b8`.
Configuration: `configs/external/comparison-audit-v1.json`.
All 6,525 regression compound/endpoint labels are scored for each variant. The
adopted baseline is a performance reference with a different feature/training
pipeline, not an isolated architecture ablation.

## Results

Lower ST-RAE and MAE are better; higher Spearman is better. Macro metrics average
the four endpoint metrics computed from pooled OOF, not fold metrics.

| Model | MA-ST-RAE | MA-MAE | MA-Spearman |
|---|---:|---:|---:|
| Adopted TabPFN 2.6 reference | 0.724704 | 0.559647 | 0.568038 |
| Matched neural control | 0.824695 | 0.609828 | 0.512470 |
| Matched neural + external auxiliaries | 0.797301 | 0.596715 | 0.517755 |

| Endpoint | Reference ST-RAE | Control ST-RAE | External ST-RAE | Control Spearman | External Spearman |
|---|---:|---:|---:|---:|---:|
| CYP1A2 | 0.816778 | 0.918743 | 0.897863 | 0.432546 | 0.437882 |
| CYP2C9 | 0.625471 | 0.719540 | 0.685803 | 0.598314 | 0.608756 |
| CYP2D6 | 0.928701 | 1.023734 | 1.013723 | 0.297603 | 0.296836 |
| CYP3A4 | 0.527865 | 0.636761 | 0.591816 | 0.721417 | 0.727545 |

All four pooled ST-RAE point estimates improve relative to the matched control,
but no endpoint beats its adopted TabPFN reference. CYP2D6 ranking does not improve;
its ST-RAE remains above the metric's evaluated-population mean reference.

An appreciable part of the benefit is in fold 3: CYP2C9 ST-RAE moves
0.744124 → 0.594369, and CYP3A4 0.741215 → 0.534943. Improvements are not equally
large across folds. In fold 3, the external model's selected epoch is 27, whereas
selection generally occurs much earlier in other folds. Identical early-stopping
rules can select different training durations; the recipe comparison includes that.

## Conditional uncertainty

Paired bootstrap draws the fixed UMAP groups with replacement, with the same group
draw applied to both models and all endpoints. It recomputes each resample's
ST-RAE denominator. There are 2,000 draws, seed 20260911. Negative deltas favour
external training.

| Endpoint | 2.5%–97.5% interval for external minus control ST-RAE |
|---|---:|
| CYP1A2 | [-0.074586, +0.014759] |
| CYP2C9 | [-0.076352, -0.004153] |
| CYP2D6 | [-0.031479, +0.008696] |
| CYP3A4 | [-0.124264, +0.000278] |
| Macro | [-0.074801, +0.000153] |

These are conditional intervals for already fitted OOF predictions. They do not
include retraining/seed uncertainty, do not prove causality for a particular source,
and are not adjusted for comparing multiple endpoints. The macro interval crosses
zero slightly. CYP2C9 has the clearest conditional signal in this run, but neither
the dataset's size nor this single comparison warrants a strong general claim.

## Run lineage and checks

| Fold | Control run | External run |
|---|---|---|
| 0 | `20260911T053158423216Z-9c04a430` | `20260911T053308971388Z-587ff71b` |
| 1 | `20260911T053556181042Z-b9378cc2` | `20260911T053658911105Z-1ea68598` |
| 2 | `20260911T053904509268Z-fa9c30a0` | `20260911T054002259749Z-29496a21` |
| 3 | `20260911T054158295368Z-f13076ec` | `20260911T054310514050Z-221b89e2` |
| 4 | `20260911T054826374242Z-0d93cdb3` | `20260911T054927294591Z-30e4d318` |

The audit verifies OOF unique keys, exact label agreement, complete coverage,
expected folds, identical matched memberships, train/validation/outer connectivity
separation, finite values, and all saved fold scores against the official scorer.
Each model's test output has exactly the expected 750 compound IDs × four
endpoints; fold-member predictions and five-fold means are retained in the audit
run with explicit compound/endpoint/split identity. They are not submission files.

The audit outputs `metrics.csv`, `macro.csv`, `paired-cluster-bootstrap.csv`,
variant-specific OOF tables, test member/ensemble tables and `audit.json`.
Registration `20260911T055252550746Z-98a54b81` appends ten model evaluations,
13,050 new OOF prediction records and the fold/pooled metrics to PostgreSQL.

Initial audit `20260911T055142656046Z-15ee49fc` failed on NumPy/Polars finite-value
interop; explicit NumPy conversion fixed the audit and the successful run
supersedes it. No models were retrained for this audit correction. Earlier retained
acquisition, curation and launch failures are listed in the external-data report.
Infrastructure checks and 39 tests, including the database integration tests, pass.

## Decision

Retain the curated external asset and the comparison as evidence; reject this
neural recipe as a replacement submission. Keep the existing TabPFN 2.6 baseline.
A next discriminating experiment is to use external training for representation
learning before a TabPFN 2.6 readout, with fold-safe fitting, or to separate external
pretraining from challenge adaptation. This has not been run here. Source-specific
ablations would be needed to attribute effects to ChEMBL versus PubChem.
No calibration, classification optimization, or new submission was performed.
