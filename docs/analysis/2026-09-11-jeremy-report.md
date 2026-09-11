# Jeremy CYP reproduction kit: source review

Reviewed public commit `d3a5549e0c98167435cbf8f5115c68a9a184ed7a` of
`jeremycheminf/openadmet_scripts`. Recorded source-fetch run:
`20260911T143603169179Z-89f6e3fd`. Source files, tree and SHA-256 provenance are
preserved in that run. The user-authorized ghq clone is at
`/home/nagaet/ghq/github.com/jeremycheminf/openadmet_scripts` and has the same HEAD.
No third-party training code was executed.

## Main finding: a different calibration target

The README reports kit OOF macro ST-RAE .614 and broader research ensemble .604.
It reports live-board .7247 with own-population placement and .5214 after
blind-moments calibration. The paragraph below the result table explicitly
attributes .5214 to the broader research ensemble. It is not independently
verified here as the exact trimmed kit's live score.

`src/cyp_submission/blind_benchmark.py` imports public constants attributed to
SuperCowPowers/briford, inferred from three affine-transformed, scored submissions:

| Endpoint | Inferred live-label mean | Inferred live-label SD |
|---|---:|---:|
| CYP1A2 | 4.412 | 1.553 |
| CYP2C9 | 4.830 | 1.101 |
| CYP2D6 | 3.107 | 1.599 |
| CYP3A4 | 4.880 | 1.272 |

These are author-reported inferred aggregate statistics, not released test labels.
The cited SuperCowPowers script contains the same numbers and the reverse-solving
code. This review checked that source but did not independently solve the moments.

`07_build_submission.py` applies, per endpoint:

    placed = mean_blind + ((pred_test - mean(pred_test)) / sd(pred_test))
                          * rho_oof * sd_blind

The rank order is unchanged for positive rho. This changes both location and
spread of final predictions. Unlike our previous OOF-derived calibration, its
location/scale target comes from live-board feedback. The README describes its
own training-population 2D6 mean around 4.8 versus inferred live mean around 3.1.
The constants refer to the currently scored subset, not necessarily the entire
blind population. This is material to what the .5214 result means.

The current upstream SuperCowPowers source also explicitly distinguishes R2-optimal
placement from ST-RAE-optimal placement: reported 2D6 ST-RAE worsened from .565 to
.694 when moving to the R2 optimum. It includes a separately board-probed 2D6
placement (mean 3.57, SD .90), while Jeremy's kit uses the simpler blind-moments
formula. Therefore the formula is a useful candidate, not a proven ST-RAE optimum.

## Modeling differences confirmed in code

| Area | Jeremy kit | Current local comparison |
|---|---|---|
| Classical model | ECFP4 + 217 RDKit2D; LGBM default L2, 500 trees, 31 leaves | FP + Mordred; L1, 300 trees, 15 leaves |
| Embeddings | Two external CYP/ADME-pretrained checkpoints: chemprop_medium and chemprop_chemeleon | Original frozen CheMeleon foundation |
| Reduction | Median imputation, StandardScaler, PCA200 | Emb-only component uses raw 2048 dimensions |
| TabPFN | 8 estimators; model version unpinned | Explicit 2.6, 8 estimators |
| Extra model | Direct 17-head Chemprop from chemprop_medium | Separate log2fc auxiliary predicts features for downstream pIC50 |
| Ensemble | Raw and affine-calibrated candidates; bagged Caruana selection | Prespecified equal pairs/triple |
| Test training | LGBM/TabPFN all labeled rows refit | Submitted TabICL is five-fold average |

The adapted `chemprop_chemeleon` is not the vanilla CheMeleon checkpoint. Its source
comments explicitly say CheMeleon initialization followed by CYP/ADME pretraining.
The README reports vanilla CheMeleon redundant with classical descriptors, agreeing
qualitatively with our observed high A/B prediction correlation. This does not
establish that every external-data training recipe will help.

The direct model has four challenge pIC50 heads, four measured log2fc heads, four
ChEMBL heads and five PubChem qHTS heads; auxiliary weights .3 and primary inverse
label-count weights. External measurements remain separate tasks. Unlike the
all-data tabular refits, final Chemprop retains one fifth of challenge data for
validation. Thus not all model families use every label for gradient updates.

Caruana implementation: 20 random candidate-subset bags, half the candidate
library per bag, initialization with the best member, up to 50 greedy additions
with replacement, ST-RAE as objective. Aggregate selection counts determine
nonnegative weights. Although README says five regression models, the submission
script lists four base sources, each with raw and calibrated variants.

## What the quoted CV number measures

Their Butina uses Morgan similarity cutoff .6 (distance threshold .4), with
cluster-size balancing and no local anchor-retention protocol. The OOF populations
and splits therefore differ from our current comparison.

PCA/scaling is fitted to all training structures before splitting, including
validation structures, although labels are not used by that transform. More
importantly, affine fits and Caruana weights are estimated on pooled OOF labels
and their final score is reported on those same labels. That is an ensemble
selection score rather than an additional held-out evaluation of the ensemble.
Our equal-blend .6736 cannot be compared as an identical evaluation to their .614.

## Practical next comparisons

1. Reuse saved OOF predictions to compare Caruana selection against equal blends,
   with meta-level held-out folds for fitting affine corrections/weights.
2. Separately inspect a blind-moment-adjusted version of current test predictions;
   save explicit constants and source provenance. Do not describe its benefit as
   a CV improvement or silently apply it to OOF labels from another population.
3. Evaluate their CYP/ADME-adapted frozen checkpoints on our matched folds, with
   train-fold StandardScaler/PCA200. This is a concrete encoder difference, not
   another run of our existing vanilla embedding model.

No candidate or submission was modified in this review.

## Sources

- https://github.com/jeremycheminf/openadmet_scripts/tree/d3a5549e0c98167435cbf8f5115c68a9a184ed7a/CYP_Challenge
- https://github.com/jeremycheminf/openadmet_scripts/blob/d3a5549e0c98167435cbf8f5115c68a9a184ed7a/CYP_Challenge/scripts/07_build_submission.py
- https://github.com/jeremycheminf/openadmet_scripts/blob/d3a5549e0c98167435cbf8f5115c68a9a184ed7a/CYP_Challenge/src/cyp_submission/tabular_cv.py
- https://github.com/SuperCowPowers/workbench/blob/main/ml_pipelines/OpenADMET/cyp/scripts/cyp_recalibrate.py (live source checked during this review)
