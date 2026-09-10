# TabPFN 3.0 comparison on both tracks

The user requested testing version 3.0 for classification and then also regression.
Use the same PCA256 + nested predicted log2fc4 + Mordred recipe from RL-0008,
with the same UMAP folds, eight estimators, seeds 42 + fold and classification
threshold 0.5. Reuse the exact 260-dimensional source matrices; train-only Mordred
filtering and median imputation are unchanged. No feature search or threshold tuning.

The Python package remains locked at 8.5.0. Select `ModelVersion.V3` explicitly and
snapshot the existing local v3 default classifier/regressor checkpoints into each
run. This compares version-specific factory configurations and weights; it is not
a claim that only tensor weights changed while internal preprocessing stayed fixed.
Effective constructor parameters are saved for each endpoint. Base weights have
run-manifest SHA-256 hashes, together with the source/config/environment inventory.

The variant is `stack_mordred_v3`, launched with
`scripts/launch_model_screen.py stack_mordred_v3`. Scoring and database validation
use `scripts/finalize_model_screen.py stack_mordred_v3`. All six endpoint models
are fitted independently in each outer fold. Candidates and submissions are not
produced. Prior v2.6 run artifacts remain unchanged.

## Coverage follow-up

The v3 classifier emitted a feature-coverage warning: with eight estimators and
200 features per member, some retained features are never sampled (e.g. 1,639
input columns exceeds 1,600). Preserve the eight-estimator fixed-budget comparison
and additionally run nine estimators on both tracks (`stack_mordred_v3_9`). This
is an explicitly separate configuration, not a replacement for the original run.
Changing estimator count can change predictions beyond feature coverage alone.

## Checkpoints and eight-estimator results

Local default checkpoint SHA-256 values:

- Regressor: `311ce18d97e9533d8585eaadafe040fbdd8070533209ed8696641dadc97a7301`
- Classifier: `d0d865d54dfbc524f5703104be90620182dca7e5fb2c16de72e9959ea18f3988`

Eight-estimator audit `20260910T130329552247Z-12eee881` verifies all 11,606 unique
OOF pairs and official scoring. Macro ST-RAE 0.730789, MAE 0.561217,
MCC 0.111502, AUROC 0.649440. This is worse than the same feature recipe with
v2.6 on all four aggregate metrics. CYP2D6 TDI stays all-negative at 0.5.
DB load: `20260910T130330644986Z-d7fe5d15`.

Eight-estimator fold runs:

- `20260910T130031103536Z-1385427e`
- `20260910T130107420355Z-af1774d5`
- `20260910T130137574127Z-b08a973c`
- `20260910T130208134152Z-72c9209d`
- `20260910T130238316435Z-b65bcb4b`

## Nine-estimator result and decision

Audit `20260910T130711842177Z-90a5a460` passed for all 11,606 OOF pairs, with
DB load `20260910T130712966180Z-79d1dc5d`. No feature-coverage warnings occurred
in any of the five nine-estimator run logs. Metrics did not improve.

| Version / estimators | Macro ST-RAE ↓ | MAE ↓ | MCC ↑ | AUROC ↑ |
|---|---:|---:|---:|---:|
| 2.6 / 8 (RL-0008) | 0.724704 | 0.559647 | 0.129339 | 0.668071 |
| 3.0 / 8 | 0.730789 | 0.561217 | 0.111502 | 0.649440 |
| 3.0 / 9 | 0.731002 | 0.561399 | 0.108931 | 0.648998 |

Nine-estimator fold runs:

- `20260910T130336415446Z-d855f196`
- `20260910T130410085602Z-3abaea33`
- `20260910T130447710493Z-e7a80e44`
- `20260910T130523774547Z-e8606e87`
- `20260910T130600596333Z-193ce1c7`

Keep v2.6 for this feature recipe on both tracks. Both v3 configurations worsen
all four aggregate metrics, and CYP2D6 TDI stays all-negative. The coverage follow-up
does not rescue the result, though it is not an isolated causal test of coverage.
This fixed-split, fixed-seed comparison does not establish that v3 is generally
inferior across feature sets or tuning settings. All ten model runs completed
without failures; saved fit states, predictions and negative evidence are retained.
No new submission or threshold optimization. Infrastructure checks and 32 tests
passed during implementation.
