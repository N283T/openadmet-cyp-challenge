# LightGBM split comparison

Pooled OOF, endpoint macro averages; ST-RAE lower is better, MCC higher is better.

| Split | Features | ST-RAE | MAE | MCC | AUROC |
|---|---|---:|---:|---:|---:|
| random | mordred | 0.7466 | 0.5766 | 0.1589 | 0.6786 |
| random | mordred_morgan | 0.7339 | 0.5700 | 0.1814 | 0.6857 |
| scaffold | mordred | 0.7506 | 0.5784 | 0.1812 | 0.6739 |
| scaffold | mordred_morgan | 0.7389 | 0.5723 | 0.1738 | 0.6859 |
| umap | mordred | 0.7658 | 0.5865 | 0.1383 | 0.6365 |
| umap | mordred_morgan | 0.7534 | 0.5802 | 0.1652 | 0.6479 |

## Mean nearest-neighbor Tanimoto

| Endpoint | Split | Mean max Tanimoto | Fraction >0.5 | Min fold n | Max fold n |
|---|---|---:|---:|---:|---:|
| CYP1A2_pIC50_direct_inhibition | random | 0.3912 | 0.1204 | 265 | 312 |
| CYP1A2_pIC50_direct_inhibition | scaffold | 0.3841 | 0.1098 | 269 | 293 |
| CYP1A2_pIC50_direct_inhibition | test_vs_all_labeled | 0.4924 | 0.5453 | 750 | 750 |
| CYP1A2_pIC50_direct_inhibition | umap | 0.3130 | 0.0099 | 228 | 379 |
| CYP2C9_pIC50_direct_inhibition | random | 0.4221 | 0.1486 | 237 | 289 |
| CYP2C9_pIC50_direct_inhibition | scaffold | 0.4238 | 0.1502 | 245 | 274 |
| CYP2C9_pIC50_direct_inhibition | test_vs_all_labeled | 0.4807 | 0.5453 | 750 | 750 |
| CYP2C9_pIC50_direct_inhibition | umap | 0.3364 | 0.0140 | 247 | 268 |
| CYP2D6_is_TDI | random | 0.3642 | 0.0949 | 283 | 318 |
| CYP2D6_is_TDI | scaffold | 0.3616 | 0.0882 | 285 | 308 |
| CYP2D6_is_TDI | test_vs_all_labeled | 0.4639 | 0.4547 | 750 | 750 |
| CYP2D6_is_TDI | umap | 0.3073 | 0.0094 | 247 | 327 |
| CYP2D6_pIC50_direct_inhibition | random | 0.3643 | 0.0951 | 282 | 316 |
| CYP2D6_pIC50_direct_inhibition | scaffold | 0.3618 | 0.0884 | 284 | 306 |
| CYP2D6_pIC50_direct_inhibition | test_vs_all_labeled | 0.4639 | 0.4547 | 750 | 750 |
| CYP2D6_pIC50_direct_inhibition | umap | 0.3072 | 0.0094 | 247 | 326 |
| CYP3A4_is_TDI | random | 0.5241 | 0.4322 | 694 | 740 |
| CYP3A4_is_TDI | scaffold | 0.4923 | 0.3792 | 705 | 729 |
| CYP3A4_is_TDI | test_vs_all_labeled | 0.5234 | 0.6227 | 750 | 750 |
| CYP3A4_is_TDI | umap | 0.3539 | 0.0360 | 655 | 761 |
| CYP3A4_pIC50_direct_inhibition | random | 0.4718 | 0.3075 | 440 | 491 |
| CYP3A4_pIC50_direct_inhibition | scaffold | 0.4595 | 0.2857 | 453 | 487 |
| CYP3A4_pIC50_direct_inhibition | test_vs_all_labeled | 0.5205 | 0.6187 | 750 | 750 |
| CYP3A4_pIC50_direct_inhibition | umap | 0.3452 | 0.0278 | 417 | 529 |

One seed, one learner, no tuning; no claim of uncertainty or leaderboard calibration.
Validation uses 80% reference labels; the test comparison uses all endpoint labels.

## Decision and limits

**Historical recommendation, superseded by [RL-0003](../research/2026-09.md#rl-0003--umap-only-and-local-research-log):**
Routine evaluation now uses UMAP only; reconsider if LB substantially diverges from CV.
The original recommendation below is retained as research history.

Use UMAP as a provisional fixed benchmark for model comparisons, following PXR,
with random as a secondary reference. Scaffold and random have similar baseline
scores; do not spend further tuning budget on choosing a split at this stage.

UMAP is substantially more structurally distant than the blinded test population:
mean maximum Tanimoto 0.307–0.354 versus test 0.464–0.523. Random is 0.364–0.524.
The test reference pool uses all endpoint labels, versus 80% in CV, so this is not
an equal-training-size comparison. In particular, UMAP should not be presented as
a calibrated proxy for leaderboard scores. No test labels were used.

Morgan improves regression macro ST-RAE by about 0.012 in all three splits.
Retain Mordred alone as an ablation and Mordred + Morgan as a stronger current
baseline. Classification gains are inconsistent across splits, and CYP2D6 TDI is
weak (UMAP + Morgan AUROC 0.536, MCC 0.094). No threshold tuning was attempted.
Next investigate predicted log2fc and frozen auxiliary embeddings with TabPFN
readout, using the fold-exclusion and inner-cross-fitting protocol in
[split-baselines.md](../split-baselines.md).

## Provenance and audit

- Preprocessed dataset: `20260909T132734119233Z-a7063ccd`.
- Mordred descriptors: `20260909T132918829974Z-bc81af58`.
- Split assignments: `20260909T135130268493Z-edba6fc0`.
- LightGBM models and OOF: `20260909T135535162532Z-eb0dd9c8` (180 models).
- Official-score audit and authoritative metrics: `20260909T140626284108Z-b5472b1e`.
- 69,636 keyed OOF rows; 216 endpoint/strategy/features/fold score groups validated.

The first baseline run's metrics Parquet inferred its schema from the first 100
regression rows and omitted classification metric columns. The JSON summary and
OOF predictions were intact. The separate audit recomputed all metrics from OOF,
verified pooled values against the JSON summary and ST-RAE against the pinned
[official scorer](https://github.com/OpenADMET/CYP-Challenge-Tutorial/blob/858ae63ce79934113bccdb7fc65467de5f7b1935/evaluation/custom_scoring_functions.py),
and saved complete metrics without changing the original run. The writer is now
fixed to infer all rows, with a mixed-track schema regression test. DB metrics
refer to the audit run; model/OOF records refer to the original baseline run.

All three split candidates remain available. No auxiliary pretraining, TabPFN fit,
blinded-test prediction, or submission was performed.

Current decisions and next experiments: [research log](../research/README.md).
Historical discussion: [Issue #1](https://github.com/N283T/openadmet-cyp-challenge/issues/1).
