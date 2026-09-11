# Regression baseline diagnosis before external-data expansion

Recorded 2026-09-11. Classification work is deferred by user decision. This is a
regression diagnosis, not a new model comparison or submission.

The baseline's pooled OOF predictions are nearly mean-calibrated, despite weak
ranking for CYP2D6 and substantial shrinkage at the label extremes. This supports
prioritizing additional training evidence over another round of readout settings.
It does **not** identify the blinded population's correct centre or prove that
calibration cannot improve leaderboard ST-RAE.

## Evidence and reproducibility

Diagnostic run: `20260911T050944039369Z-36d00d9c` (regression, completed).
Source: `scripts/regression_diagnostics.py`; configuration:
`configs/analysis/regression-diagnostics.json`. Its manifest declares the processed
measurements, audited OOF/metrics, test ensemble/member predictions, fixed UMAP
assignments, and official scoring implementation. The run snapshots all of them.

- Dataset: `20260909T132734119233Z-a7063ccd`.
- Audited OOF: `20260910T123856868493Z-243a1e23`.
- Test inference: `20260910T132141471698Z-343eb5ec`.
- Split: `20260909T135130268493Z-edba6fc0:umap`.
- Existing LB evidence is documented in [baseline submission](../submission-baseline.md).

Audits passed: unique compound/endpoint keys; all 6,525 regression labels covered;
labels exactly agree with OOF truth; recorded folds match the fixed assignments;
no labeled/test compound overlap; 750 test compounds per endpoint; five predictions
per test compound; member means reproduce the saved ensemble; finite values and
ordered intervals; all pooled and fold ST-RAE, MAE and Spearman reproduce the
previous audit. ST-RAE is evaluated with the declared official scoring source.

Outputs under `artifacts/runs/regression/20260911T050944039369Z-36d00d9c/outputs/`:
`summary.json`, `metrics.csv`, `distributions.csv`, `residual_bins.csv`,
`regression-oof-diagnostics.parquet`, `test_disagreement.csv`,
`test_member_correlations.csv`, keyed per-endpoint test-disagreement Parquets,
and `regression-diagnostics.png`. Residual sign is prediction minus observation;
SD uses the population convention (ddof=0). Prediction tables retain compound,
endpoint and split/fold identity.

## OOF behaviour

| Endpoint | Labels | ST-RAE | Spearman | Mean residual | SD(pred)/SD(label) | Fold ST-RAE range |
|---|---:|---:|---:|---:|---:|---:|
| CYP1A2 | 1,412 | 0.816778 | 0.501103 | +0.010954 | 0.510 | 0.768–0.897 |
| CYP2C9 | 1,285 | 0.625471 | 0.648118 | +0.002214 | 0.613 | 0.492–0.699 |
| CYP2D6 | 1,493 | 0.928701 | 0.360647 | -0.002295 | 0.365 | 0.871–1.011 |
| CYP3A4 | 2,335 | 0.527865 | 0.762283 | -0.006795 | 0.754 | 0.479–0.605 |

Across all endpoint-specific predicted quintiles, mean residuals lie between
-0.065 and +0.071 pIC50. There is no large pooled offset or obvious monotonic
mean-residual trend to correct. This is descriptive, without uncertainty intervals;
it does not establish perfect calibration, especially for an absolute-error metric.

CYP2D6 is weak throughout the folds: Spearman ranges from 0.311 to 0.443.
Its fold-0 mean residual is -0.211, while fold 2 is +0.130, so pooled near-zero
bias hides some chemistry-dependent variation. Fold 3 ST-RAE exceeds 1.0.
The score's reference mean is computed on that evaluated fold, not a deployable
constant fitted on its training set. Fold variability is not seed variability or
an independent estimate of statistical uncertainty.

All four endpoints have zero OOF predictions at or above pIC50 6, although their
labels contain such compounds. In CYP2D6, 129 labels below 4 have mean truth 2.679
and mean prediction 4.518; 97 labels in [6,7) have mean truth 6.423 and mean
prediction 4.976. These show a failure to separate extremes. Grouping by the true
label cannot provide a usable test-time calibration rule. A narrow prediction
range alone does not justify widening: weak predictors should shrink toward the
conditional mean, and widening can increase errors on incorrectly ranked compounds.

## Soft-threshold error

| Endpoint | Mean absolute error | Mean error outside label CI | Predictions inside label CI |
|---|---:|---:|---:|
| CYP1A2 | 0.652 | 0.401 | 22.9% |
| CYP2C9 | 0.451 | 0.196 | 43.6% |
| CYP2D6 | 0.599 | 0.389 | 19.0% |
| CYP3A4 | 0.537 | 0.277 | 38.2% |

For CYP2D6 labels below 4, MAE is 1.839 but mean error outside the interval is
0.707. For labels in [6,7), the corresponding values are 1.447 and 1.334.
Thus large errors in opposite activity regions receive different effective
penalties. Optimizing mean/RMSE calibration is not equivalent to optimizing ST-RAE.
CI coverage here describes whether a point prediction lies in the measurement's
credible interval; it is not coverage of a model uncertainty interval.

## Blinded prediction distributions

| Endpoint | Label mean ± SD | OOF mean ± SD | Test ensemble mean ± SD |
|---|---:|---:|---:|
| CYP1A2 | 4.955 ± 1.030 | 4.966 ± 0.525 | 4.947 ± 0.601 |
| CYP2C9 | 4.581 ± 0.782 | 4.583 ± 0.479 | 4.745 ± 0.590 |
| CYP2D6 | 4.784 ± 0.916 | 4.782 ± 0.334 | 4.593 ± 0.312 |
| CYP3A4 | 4.096 ± 1.093 | 4.089 ± 0.824 | 4.551 ± 0.838 |

CYP2D6's test ensemble central 90% is 4.104–5.163; only 10/750 predictions are
below 4. Individual member SDs are 0.273–0.360, so the concentration is present
before averaging. Member mean predictions range from 4.515 to 4.640. Averaging
contributes to compression, but is not its sole source.

CYP3A4 test predictions shift upward relative to OOF; CYP2C9 also shifts upward,
while CYP2D6 shifts downward. These are shifts in **model predictions**, not
observed test labels. They cannot show whether the shifts are sufficient or correct.

The September 10 LB snapshot shows CYP2D6 MAE 1.724 versus OOF 0.599, with similar
Spearman (0.374 versus 0.361). Population/calibration mismatch remains plausible,
but its direction and magnitude are not identified by MAE and rank correlation.
No hidden-label centre is inferred or borrowed from the community article.

## Existing ChEMBL source and next experiment

The prior PXR project has a working local PostgreSQL ChEMBL 36 database:
`chembl_36`, Unix socket `/tmp`, port 5433. The socket location is not the location
of experiment outputs. The prior archive is under
`/home/nagaet/pxr-iduction-challenge/data/chembl/chembl_36/chembl_36_postgresql`.
The database reports ChEMBL_36, created 2025-07-28; the community article used 37.

Read-only, repeatable-read inventory run:
`20260911T051139350452Z-8a9e87f9` (shared, completed).
Query/config: `scripts/chembl36-inventory.sql`,
`configs/external/chembl36-inventory.json`. Result snapshots include measurement
strata, human CYP target metadata, executed SQL and database version/provenance.
The first attempt, `20260911T051112042187Z-ad931847`, failed while serializing a
metadata timestamp after the query succeeded. It remains recorded and is superseded
by the successful run; no source database changes were made.

| Endpoint | Raw IC50 activity rows | IC50 rows with pChEMBL | Raw AC50 activity rows |
|---|---:|---:|---:|
| CYP1A2 | 5,684 | 1,672 | 17,143 |
| CYP2C9 | 7,344 | 3,027 | 17,143 |
| CYP2D6 | 8,493 | 3,279 | 17,143 |
| CYP3A4 | 13,587 | 6,494 | 17,143 |

These are uncurated measurement rows, not unique molecules or additional usable
labels. They include differing relations, assay types and quality flags.
Missing pChEMBL does not itself mean inactivity. Equal AC50 row counts across
isoforms warrant checking assay/document provenance before combining this source
with a separately downloaded qHTS panel; overlap is not yet verified.

Next: extract source-preserving records with assay/document identity, units,
relations, quality flags and standardized structures; verify possible qHTS overlap;
exclude held-out molecular identities from external training by fold. Keep external
readouts separate from challenge pIC50 labels. Then compare challenge-only and
external-augmented training under the same regression UMAP protocol. Calibration
remains a separate later comparison, rather than being coupled to the data change.
Classification tuning is deferred until regression stalls or priorities change.
No new model training, external-data integration, calibration or submission has
been performed in these diagnostic/inventory runs.
