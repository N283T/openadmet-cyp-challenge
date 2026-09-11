# Separate-feature ensemble components

The user proposed an ensemble of three separate regression models. Earlier
embedding-only OOF predictions used different folds and a fine-tuned encoder,
so they are not pooled with current folds. This experiment trains all three
components on the current parent-retained series folds and adopted auxiliary
configuration.

- Morgan radius-2, 2,048-bit fingerprints (default non-chiral Morgan) plus Mordred
  descriptors -> LightGBM L1, 300 trees, 15 leaves, minimum child size 20,
  learning rate .05, L2 1. Missing/constant Mordred columns and median imputation
  are determined from each endpoint's training rows only.
- Frozen CheMeleon 2,048-dimensional message-passing embeddings -> TabPFN 2.6.
  No external PCA or predicted log2fc is added to this component.
- Four adopted nested predicted-log2fc values -> TabPFN 2.6. Each outer/recipient
  fold pair is excluded from the corresponding auxiliary training. Outer-fold
  prediction features average four eligible auxiliary models.

Both TabPFN components use eight estimators and seed 42 + fold. The prediction
response is direct-inhibition pIC50. All components retain anchors in training,
exclude anchors from validation, and use identical validation molecules.
Classification, test prediction and submission are outside this comparison.

The audit reports per-endpoint prediction Pearson/Spearman correlations and
Pearson correlations of residuals (prediction minus observed pIC50), on pooled
OOF and within individual folds. It also evaluates each component, all three
pairwise equal means, and the three-way equal mean with the official ST-RAE,
MAE and Spearman scores. No ensemble weights are optimized on these OOF labels.
Current full-feature TabICL/TabPFN OOF predictions provide matched controls.
Neighbor-only results are retained separately.

Initial run `20260911T142309831305Z-0b7fb03a` completed model fitting but failed
when concatenating float32 and float64 prediction columns. The failed run and
saved fits remain intact. The output prediction column now explicitly uses
float64 for every component; the complete fold is rerun.

## Completed runs and results

Successful folds 0–4:
`20260911T142441093420Z-91e36a41`, `20260911T142530231264Z-a5b5b04c`,
`20260911T142619356692Z-5c334eef`, `20260911T142708466901Z-8b7b5050`,
`20260911T142757006049Z-db34f643`.
All 60 endpoint/component/fold fits and their model assets were saved.
Audit `20260911T142905405713Z-d0f742da` checks 6,384 unique validation keys per
component and exact label/bound/fold agreement with current controls. It saves
aligned OOF predictions, per-fold and pooled correlations, scores, and macro
summaries. Scores are also registered in PostgreSQL.

Let A = FP+Mordred LightGBM, B = frozen embedding TabPFN, C = pred_log2fc TabPFN.
Each correlation below is pooled OOF Pearson, followed by residual Pearson in
parentheses. Residuals are prediction minus raw pIC50.

| Endpoint | A vs B | A vs C | B vs C |
|---|---:|---:|---:|
| 1A2 | .884 (.965) | .772 (.918) | .825 (.931) |
| 2C9 | .867 (.929) | .795 (.857) | .844 (.884) |
| 2D6 | .801 (.971) | .701 (.961) | .740 (.958) |
| 3A4 | .916 (.904) | .851 (.809) | .898 (.866) |

A vs C has the lowest prediction correlation for all four endpoints. Residual
correlation remains high, particularly for 2D6; these models still share much of
their error. Correlation alone does not select an ensemble; direct blend scores
below demonstrate the actual improvement on these folds.

| Method | All OOF macro ST-RAE | Neighbor macro ST-RAE |
|---|---:|---:|
| A | .728820 | .858036 |
| B | .721992 | .895810 |
| C | .698059 | .940962 |
| (A+B)/2 | .704776 | .853509 |
| (A+C)/2 | **.673634** | .847109 |
| (B+C)/2 | .677522 | .876974 |
| (A+B+C)/3 | .674589 | .843485 |
| Current full TabPFN | .688696 | .854329 |
| Current full TabICL | .695341 | **.834264** |

| Endpoint | (A+C)/2 | Three-way mean | Current full TabPFN | Current full TabICL |
|---|---:|---:|---:|---:|
| 1A2 | .757879 | .756238 | .775548 | .767708 |
| 2C9 | .561990 | .567415 | .570962 | .585560 |
| 2D6 | .877768 | .874623 | .900283 | .920769 |
| 3A4 | .496899 | .500080 | .507991 | .507326 |

Both A+C and the three-way mean improve all four endpoint ST-RAEs relative to
both current foundation controls on all OOF. A+C macro MAE is .533516 and
macro Spearman .591016, versus current full TabPFN .536118 and .589681.
The three-way mean has MAE .533836 and Spearman .593827.

Decision: retain A+C equal mean and the three-way equal mean as simple ensemble
candidates. Their all-OOF performance supports separating the feature families
rather than relying only on feature concatenation. The neighbor subset does not
beat the current full TabICL ST-RAE. These are CV findings, with no new LB result;
no weights were optimized and no automatic current-model pointer replacement or
submission was made.
