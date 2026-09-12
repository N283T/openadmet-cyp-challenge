# CYP2D6 specialist scope and initialization comparison

This adapts the specialist-encoder idea in the updated
[SuperCowPowers working notes](https://supercowpowers.github.io/workbench/blogs/cyp_challenge/).
The source describes encoders with different target sharing. The present trial
isolates challenge-head scope after a common external initialization; it does not
reproduce the source's joint six-head external-data training or placement correction.
External ChEMBL/qHTS information enters through the same Jeremy medium pretrained
encoder as existing experiments. A separate external-efficacy probe is recorded
in the parallel experiment branch.

## Fixed design

Five arms use the same 600-dimensional graph architecture, fixed parent folds,
and CYP2D6-only inner ST-RAE checkpoint selection:

| Arm | Challenge outputs | Initialization | Encoder/head LR |
|---|---|---|---|
| Full | Direct, TDI, Emax, log2fc for all four CYPs (16) | External pretrained | 1e-5 / 1e-4 |
| 2D6 auxiliary | 2D6 direct, TDI, Emax, log2fc (4) | External pretrained | 1e-5 / 1e-4 |
| 2D6 direct | 2D6 direct pIC50 only (1) | External pretrained | 1e-5 / 1e-4 |
| Fresh, matched LR | 2D6 direct only (1) | Random | 1e-5 / 1e-4 |
| Fresh, scratch LR | 2D6 direct only (1) | Random | 1e-3 / 1e-3 |

Direct outputs use the confidence-interval loss, weight 2. TDI, Emax, and log2fc
remain separate MSE outputs, with weights 0.5, 0.25, and 1, respectively. Missing
labels are masked, each head is normalized using inner training only, and native
signed Emax values are retained. Actual labels come from canonical measurements
and the audited observation target artifact `271e3020`; unassigned compounds
are excluded from auxiliary training. The common pretrained bridge is `7604ffeb`.

All arms use maximum 160 epochs, patience 20, batch 64, dropout 0.1, weight decay
1e-4, seed 42 plus fold, and a grouped 10% inner validation. The inner selection
criterion is CYP2D6 in every arm, including the all-CYP control. Every label of
an outer-held-out compound is excluded, and anchor-containing groups remain in
training. Changing target scope intentionally changes eligible training compounds
and shared outputs; it is not a fixed-label-count ablation.

The exact-LR fresh arm isolates initialization, while the additional scratch-rate
arm checks whether the conservative transfer learning rate unfairly handicaps
random initialization. Neither is a hyperparameter search. A previous direct-only
2D6 trial already existed; its repetition here supplies a contemporaneous control,
not a claim that direct-only modeling is new.

Each selected encoder produces all-compound 600-dimensional fold-specific
embeddings. The downstream CYP2D6 readout uses fixed LGBM (L1, 300 trees, 15 leaves,
minimum 20 samples per leaf, learning rate 0.05, L2 1) plus the new nested primary
features `d272774c`. Neural direct predictions and readout predictions are both
evaluated, without calibration or submission.

## Comparison reference

Fixed 25% and 50% additions use the latest global OOF `84c28266`, separately for
`caruana_weight_crossfit` and `caruana_full_fit`. These are different estimates;
the full-fit ensemble's score has already been used to choose weights and is not
an untouched CV benchmark. No selection-fold labels are used during base fitting.

## Results

Lower CYP2D6 ST-RAE is better. The last column is a fixed 25% addition of the
new LGBM prediction to the latest global weight-crossfit prediction (75%).

| Arm | Neural direct | LGBM readout | 25% LGBM addition |
|---|---:|---:|---:|
| Existing global reference | — | — | 0.825436 |
| Full 16 heads | 0.870812 | 0.838186 | 0.824156 |
| 2D6 auxiliary 4 heads | 0.867278 | 0.844236 | 0.823470 |
| 2D6 direct 1 head | 0.874202 | 0.838988 | 0.824223 |
| Fresh, matched LR | 0.916005 | 0.864242 | 0.824937 |
| Fresh, scratch LR | 0.908083 | 0.868519 | 0.823601 |

The 2D6-only auxiliary model improves its direct neural prediction slightly
relative to the full-head control, but its LGBM readout is worse alone. Its small
ensemble improvement is consistent with diversity rather than stronger standalone
prediction. Increasing additions to 50% gives 0.826203 for its LGBM, worse than the
existing reference. Every neural direct prediction worsens the reference at both
25% and 50% weights. Thus these trials do not show that removing cross-isoform
sharing or discarding external initialization solves the 2D6 problem.

Random initialization is weaker under both tested learning-rate settings. Fresh
models are more different from the ensemble (LGBM prediction Pearson 0.841–0.857,
versus 0.921–0.936 for transferred models), so even a weaker fresh model can make
a small contribution at low weight. The best point changes remain approximately
0.002 ST-RAE, not a large performance increase. Selection and repeated CV use
remain development comparisons, not independent confirmation of such small gains.

All ten OOF artifacts are retained for the subsequent expanded global ensemble
selection. No new production candidate, calibration, or submission is created here.

## Runs and validation

| Arm | Fit run | LGBM run |
|---|---|---|
| Full 16 heads | `20260912T101144735437Z-661a1844` | `20260912T101614614291Z-00de4aad` |
| 2D6 auxiliary | `20260912T101145931305Z-90443092` | `20260912T101620627804Z-fa092585` |
| 2D6 direct | `20260912T101147108878Z-fc9a4ac0` | `20260912T101418694424Z-390ed3c8` |
| Fresh, matched LR | `20260912T101218879142Z-b1286445` | `20260912T101606024008Z-70d8bdfa` |
| Fresh, scratch LR | `20260912T101246158987Z-b36fd27b` | `20260912T101411999351Z-601b7cb4` |

Independent comparison run: `20260912T101756235780Z-de80501c`.
It checks all ten 1,457-row CYP2D6 OOF sets against the reference by compound,
endpoint, fold, split, truth, and confidence bounds. All 25 membership files pass
outer-compound exclusion and inner-group disjointness checks. Every checkpoint
history's selection score exactly corresponds to the 2D6 inner endpoint score.
All direct predictions, fold embeddings, best/last checkpoints, histories,
memberships, and fitted LGBMs are retained. Comparison outputs include per-fold
official metrics, all prespecified blends, prediction/residual correlations, and
the verification record.

The first audit `20260912T101725660350Z-8ff89c3e` failed because the prior global
aggregate omits `split_id`. The corrected audit declares the canonical assignments,
verifies every aggregate row's fold and nonanchor membership, then attaches that
verified split ID for comparison. The original failed run is retained.

## Prespecified seed and scope means

After the initial comparison, the 2D6 auxiliary recipe is repeated with model
seeds 17 and 89, without parameter tuning. Inner split seed stays at 42, so
checkpoint validation compounds are identical across the three runs. Torch
initialization/dropout and training minibatch order vary; outer assignments and
readout random seed remain unchanged. The seed summary verifies every membership
table against seed 42 before interpreting the repetition as model-seed variation.

Two means are specified before seeing the new scores, separately for direct
neural outputs and LGBM readouts: (1) the equal mean of the 42/17/89 auxiliary
specialists; (2) the equal mean of the full-head, seed-42 auxiliary, direct-only,
and fresh scratch-rate scopes. No members are removed using their score. Both
standalone means and their fixed 25%/50% additions to each global reference are
recorded, with all constituents and transformations declared as run inputs.

| Model or fixed mean | Neural direct | LGBM alone | 25% LGBM addition | 50% LGBM addition |
|---|---:|---:|---:|---:|
| Auxiliary seed 42 | 0.867278 | 0.844236 | 0.823470 | 0.826203 |
| Auxiliary seed 17 | 0.870277 | 0.841111 | 0.822970 | 0.824695 |
| Auxiliary seed 89 | 0.873129 | 0.840989 | 0.823919 | 0.825872 |
| Three-seed mean | 0.869418 | 0.836105 | 0.823102 | 0.823991 |
| Four-scope mean | 0.856016 | 0.828992 | 0.822785 | 0.822232 |

The repeated seeds show similar performance and all three LGBM additions improve
the same reference slightly. The four-scope mean is more useful than the seed mean
here: it improves standalone direct and LGBM scores beyond its constituents and
reaches 0.822232 with the prespecified 50% LGBM addition, versus global 0.825436.
This supports modest ensemble complementarity, while leaving most of the 2D6
error unresolved. The direct neural means still worsen the global model at both
25% and 50% additions.

Seed 17 fit/readout: `20260912T101944672747Z-d726bb6d` /
`20260912T102337278051Z-f89acb5a`.
Seed 89 fit/readout: `20260912T101945866326Z-bf0ef563` /
`20260912T102341169498Z-491ace93`.
Successful seed/mean audit: `20260912T102357053293Z-e9975cb3`.
All 15 membership tables exactly match seed 42, and all seed histories select on
the 2D6 inner score. The output directory stores four mean OOF files named
`d6_seed_mean_direct-oof.parquet`, `d6_seed_mean_lgbm-oof.parquet`,
`d6_scope_mean_direct-oof.parquet`, and `d6_scope_mean_lgbm-oof.parquet`, plus the
full metrics, fixed global blends, and audit. They remain candidates for the
root experiment's global reselection rather than a separate submission.
