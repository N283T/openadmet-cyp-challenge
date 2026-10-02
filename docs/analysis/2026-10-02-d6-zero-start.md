# CYP2D6 TDI: independent restart

Subsequent status: the user resumed the investigation to connect experiments to
prior failures and carry them through to a concrete candidate. See the
[failure-linked probe](2026-10-02-d6-failure-linked-probe.md) for the matched
historical refits, unresolved causes and prepared 77-call D6 candidate. The earlier
closure below records the disposition at that time, not the current direction.

The user explicitly rejected preserving the historical D6 ensemble or its positive
call count as design anchors. This decision supersedes the October 1 recommendation
to use the historical D6 vector as the modeling control. Retained historical runs
remain evidence of failed assumptions; their fitted assets and prediction vectors
are excluded from the new modeling inputs. No rollback or new upload is authorized.

## Foundational question

Can molecular structure predict the assay's paired potency change, and can that
signal support a transferable binary decision? A gain in AP alone is insufficient
to answer this question. Neither the historical Live winner nor a fixed positive
budget establishes the right representation, validation population, or decision rule.

The canonical source has 6,146 development molecules and 750 blind molecules.
D6 has 1,497 supplied labels, including 324 positives. Of these, 1,364 have direct
pIC50 above 4 (273 positive), 51 are inferred positives with weak direct inhibition
but TDI pIC50 above 4.301, 60 are conventionally assigned inactive negatives,
and 22 fall outside these paired-measurement regimes. Missing labels are never
treated as negatives. Supplied booleans remain the scoring truth.

The official label rule is documented in the retained
[challenge snapshot](../../references/official/cyp/challenge/20260907T135802Z/snapshot.md).
The 750 blind molecules were selected through analog expansion around hits for
other enzymes, so random OOF success cannot establish Live transfer. The earlier
[acquisition analysis](2026-09-24-d6-analog-expansion-validation.md) documents this
selection problem; it is not a recipe to preserve the old calls.

## First independent experiment

Configuration: [d6-zero-start-v1.json](../../configs/classification/d6-zero-start-v1.json).
Launcher: [launch_cyp_d6_zero_start.py](../../scripts/launch_cyp_d6_zero_start.py).
Recorded run: `classification/20261002T004016766672Z-e2d420d4`.

Only canonical compounds, source measurements, and aliases are declared inputs.
Fresh Morgan fingerprints and RDKit 2D descriptors are generated from standardized
SMILES. No prior model, prediction, threshold, selected member, or positive count
is read. All fitting starts from untrained estimators.

Three independent families are compared:

- Molecular fingerprint logistic classification.
- Fingerprint plus RDKit descriptor classification with LightGBM.
- Paired-measurement learning: predict TDI pIC50 and the direct-to-TDI shift;
  score the smaller margin above their respective 4.301 and 0.301 boundaries.

Four outer random folds and four outer scaffold-group folds each contain three
inner folds. All training measurements, including auxiliary regression targets,
exclude outer held-out molecules and inner held-out molecules at the appropriate
fit. The operating threshold maximizes 1,000-bootstrap mean MCC on inner OOF
labels. There is no positive-count constraint, score-rank mixture, or prior rate.
Outer calls use those inner-selected thresholds. AP, AUROC, ordinary MCC, and
bootstrap mean MCC are recorded separately. Pooled OOF scores select a final
threshold for each family's independent blind predictions; that threshold fitting
is not counted as an unbiased performance estimate.

Scaffold groups exclude the same core scaffold across training and validation.
Acyclic molecules use individual groups, so this does not remove all close analogs.
The recorded bootstrap uses NumPy Generator draws; it follows empirical resampling
and MCC arithmetic but does not claim bit-identical official bootstrap samples.
The vectorized implementation was checked against direct sklearn MCC, including
degenerate all-positive and all-negative predictions.

The paired-measurement model is a diagnostic starting point: below pIC50 4,
point-estimate shifts can be censored or noisy. A failure must be inspected by
measurement regime and interval uncertainty before interpreting it as evidence
against learning the mechanism. Likewise, success on these folds does not resolve
the unknown scored-label mask or the blind acquisition distribution.

The initial launch `20261002T003916853752Z-2cad05be` stopped before any fit due
to an unbound variable in the random-fold constructor. Its reason is recorded;
the source audit is retained, the branch is fixed, and both split constructors
were preflight-checked before the replacement launch.

## Results

The first run completed. The retained [diagnosis report](../../artifacts/tracked/cyp-d6-zero-start-20261002T004315539799Z-6f0f98c2/report.md)
includes thresholds, regimes, and interval summaries.

| Fresh family | Random AP | Random nested MCC | Scaffold AP | Scaffold nested MCC |
| --- | ---: | ---: | ---: | ---: |
| Morgan logistic | 0.3136 | 0.1159 | 0.2869 | 0.0861 |
| Morgan/RDKit LightGBM | 0.2612 | 0.0672 | 0.2707 | 0.0636 |
| Paired measurement margin | 0.2521 | 0.0281 | 0.2647 | 0.0495 |

MCC is the bootstrap mean over outer calls produced with inner-selected thresholds,
not the maximum after fitting a threshold on pooled OOF labels. Point-estimate
paired-measurement learning is a retained negative result. Fresh molecular features
contain a weak ranking signal; regression alone is not a demonstrated TDI solution.

Among direct-active labelled molecules, conservative combinations of marginal
potency intervals put 83 entirely above the shift boundary, 378 entirely below,
and 903 across it. These are not joint posterior probabilities and do not justify
discarding or changing supplied labels. Blind call counts range from 61 to 592
and change substantially between validation strategies; no count is a target.

## Confidence and positive-regime experiment

Run `classification/20261002T004536560534Z-21c1509c` fits three more families from
the same source tables, without reading first-pilot predictions or weights.
[Configuration](../../configs/classification/d6-zero-start-confidence-v2.json).

Training examples with intervals supporting their label receive weight 1;
ambiguous or unpaired examples receive weight 0.25. Weights are normalized to mean
1 within each fitting cohort to preserve the regularization scale. The three-class
model distinguishes negatives, ordinary positives, and inferred positives, then
sums the two positive molecular probabilities. Validation uses every supplied label
with equal weight. All fitting boundaries also apply to weights and subtype targets.

| Fresh family | Random AP | Random nested MCC | Scaffold AP | Scaffold nested MCC |
| --- | ---: | ---: | ---: | ---: |
| Confidence-weighted logistic | 0.3115 | 0.0935 | 0.2941 | 0.1124 |
| Positive-regime multiclass | 0.3119 | 0.1161 | 0.2804 | 0.0560 |
| Confidence-weighted regime | 0.3092 | 0.0962 | 0.2833 | 0.0911 |

Confidence weighting modestly improves the scaffold result while lowering the random
result. Subtype separation alone does not improve structural transfer. These are
sequential hypotheses on the same development splits, not independently confirmed
improvements. They establish neither a submission choice nor a winning call count.
Models and diagnostic predictions are retained under
`artifacts/tracked/cyp-d6-zero-start-20261002T004536560534Z-21c1509c/`.

## Fresh cross-enzyme representation

Run `classification/20261002T005041934540Z-c62b6e78` trains a randomly initialized
shared molecular network on four enzymes' supplied TDI labels and paired continuous
assays. [Configuration](../../configs/classification/d6-zero-start-multiassay-v3.json).
This representation is learned from the source assay matrix. Each fitting boundary
excludes every assay measurement of held-out molecules. Missing targets are masked.
Continuous potencies are floored at 4 to avoid arbitrary sub-assay point estimates;
their difference is an activity-change auxiliary, not a measured true shift below
the reliable range. The official D6 label is the main supervised task. Architecture,
loss weights, optimizer, dropout, and 100 epochs are fixed before fitting.
Thresholds remain nested and unconstrained.

The multiassay model completes with random AP/MCC .3382/.1282 and scaffold
.3243/.1123. D6-only ablation `classification/20261002T005413496423Z-f92f58e6`
uses the same randomly initialized architecture and training budget, with auxiliary
losses exactly zero. It gives random .3383/.1514 and scaffold .3174/.1276.
Paired development intervals do not establish a gain from cross-enzyme supervision.

## Actual molecular-neighbor separation

The 1,497 D6 labels span 1,450 exact Murcko groups. That split does not guarantee
exclusion of near neighbors with different cores. Run
`classification/20261002T005927124518Z-ce5df58c` builds connected components of
the development-only Morgan graph at Tanimoto >= 0.6: 5,204 groups over 6,146
molecules, 2,077 edges, largest group 21, and 1,473 labelled groups. At both inner
and outer boundaries, no training-validation pair has Tanimoto >= 0.6. Blind
molecules do not enter clustering. Close labelled analog supervision is sparse;
grouping alone cannot create additional signal.

| Fresh family | Neighbor-excluded AP | Neighbor-excluded nested MCC |
| --- | ---: | ---: |
| Morgan logistic | 0.2907 | 0.0482 |
| D6-only network | 0.3038 | 0.1184 |
| Multiassay network | 0.3016 | 0.0881 |

Recorded comparison `classification/20261002T010208384815Z-78bd6216` uses paired
chemical-group bootstrap. D6-only network minus new logistic MCC is +.0690,
conditional 95% interval [.0174, .1239]; the AP difference interval includes zero.
Multiassay-versus-D6-only differences remain unresolved. These intervals omit
retraining and adaptive experiment selection effects. They do not establish Live
performance. See the [consolidated report](../../artifacts/tracked/cyp-d6-zero-start-20261002T010208384815Z-78bd6216/report.md).

## Experimental recommendation before the user stopped the investigation

Continue independent molecular-network development. Calibration and decision-rule
transfer remain unresolved. There is no fixed positive budget or historical-score
blend. New D6 network blind counts are 203, 94, and 227 depending on which split
supplies the final threshold; these are diagnostics, not prescribed counts.
Their sensitivity needs investigation rather than selecting a count because it
resembles a historical submission. Measurement-interval overlap is a plausible
difficulty, not a demonstrated sole cause. Negative hypotheses are retained.
No submission payload was prepared or uploaded.

## Final disposition

The user ended the investigation and requested that the changes be summarized and
pushed. Creating new molecular models was a premature interpretation of the request
to rethink the D6 problem from first principles. These trials do not establish the
cause of AP-to-Live divergence, solve decision-rule transfer, or demonstrate improvement
over the historical submitted model. The recommendation above is superseded; no
new model is adopted and no further experiment is running. Preserve the source,
configurations, fitted assets, keyed predictions and negative results as exploratory
evidence. No rollback or new submission was made during the restart.
