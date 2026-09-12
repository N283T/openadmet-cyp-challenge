# Regression OOF goal result

The final fixed five-component ensemble achieves **official macro ST-RAE
0.5964141092282933** on the retained parent CV's **6,384 OOF labels**. Its
neighbor subset score is **0.6872080091480266** on 299 labels. This reaches the
requested overall OOF range of 0.5–0.6; individual targets remain uneven.

| CYP | OOF labels | Overall ST-RAE | Neighbor ST-RAE |
|---|---:|---:|---:|
| 1A2 | 1,372 | 0.683335 | 0.726437 |
| 2C9 | 1,255 | 0.463781 | 0.777974 |
| 2D6 | 1,457 | 0.834741 | 0.832878 |
| 3A4 | 2,300 | 0.403800 | 0.411543 |

## Exact saved result and independent reconstruction

The first goal-reaching five-way mean was **0.5976166863970354** overall /
**0.6931179877186007** neighbors:

- Prediction: `20260912T072501360501Z-f8c3a807`.
- Official audit: `20260912T072502729280Z-d8928353`.
- Independent package: `20260912T073002574503Z-8be77a87`.

The completed LightGBM tuning control replaces only the fine-tuned-medium LGBM
component and improves the mean to **0.5964141092282933**:

- Prediction: `20260912T072935321335Z-1e8d2643`.
- Official audit: `20260912T072935933157Z-6e1a3bce`.
- Independent package: **`20260912T073049029140Z-0c26dd06`**.
- Exact ensemble config: `configs/external-ensemble/caruana_medium_scores_tuned_ft_cheme_v1.json`.

Each independent package rereads all five component predictions, verifies their
SHA256 hashes against the original blend manifest, reconstructs the equal mean,
and checks all molecule/endpoint/fold keys, original labels, confidence bounds and
neighbor flags against original measurements and parent assignments. Every one of
the 6,384 rows matches. Reconstruction differs by at most **8.881784197001252e-16**
from the saved predictions; it is not bitwise equal, but passes absolute tolerance
1e-12. The official metric is recomputed by endpoint, fold and subset. The final
independent macro is 0.5964141092282932, a floating-point rounding difference of
approximately 1e-16 from the original audit.

Package outputs include `README.md`, `verified-oof.parquet`, `metrics.csv`,
`macro.csv`, `membership-verification.json` and `package-manifest.json`. The
manifest contains exact weights, input hashes, original paths, run IDs and saved
component source/config/intent/result references.

## Final composition

| Component | Run | Weight |
|---|---|---:|
| Original nested Caruana ensemble | `96fce32b` | 0.2 |
| Jeremy medium frozen embedding → TabPFN | `2f0fbfbd` | 0.2 |
| Jeremy upstream FFN scores + nested log2fc → TabPFN | `07f2749a` | 0.2 |
| Jeremy medium fine-tuned embedding + nested log2fc → tuned LGBM | `e6b6cbcb` | 0.2 |
| Jeremy CheMeleon PCA256 + nested log2fc → TabICL | `73a6fb7e` | 0.2 |

The final weights were fixed equally; they were not fitted on OOF predictions.
The first goal-reaching result used `730c8bea` instead of `e6b6cbcb` in the fourth
position.

## Training provenance and interpretation

Public external pretrained weights were authorized by the user. Jeremy's public
upstream medium and CheMeleon checkpoints contain **47 pretraining output heads**;
they are distinct from his **17-head challenge fine-tune**. The public source does
not provide a training molecule manifest or named 47-task schema. This is an
unresolved provenance detail, **not evidence of leakage and not a restriction on
using the authorized weights**. The recorded provenance inspection `39e5a986`,
bridge `7604ffeb`, frozen export `c060d197` and CheMeleon export `299439e8` are
included as source references.

The independent package rereads actual five-fold fine-tuning membership from
`8d434eee`, confirming that all non-anchor outer-fold compounds are marked
`outer_holdout` and cannot enter encoder training or its inner validation. The
nested log2fc preparation `bac5fb3a` supplies the recorded ten pairwise-exclusion
membership audits. Original Caruana source manifests reference its saved grouped
inner memberships, fitted weights and outer-fit models.

The tuned LGBM performs 24 trials per endpoint/fold on grouped inner holdouts,
then refits on outer training data. Its encoder had already seen the readout's
inner-holdout labels during encoder training, so inner HPO scores can be
optimistic; this is not fully nested encoder-plus-readout selection. Outer labels
remain excluded from encoder fitting, checkpoint selection and readout HPO.
See `2026-09-12-jeremy-ft-logfc-rustuna.md` for the tuning audit.

These are development OOF results after multiple comparisons, not untouched test
or leaderboard results. The goal achieved here is the requested overall OOF
range. This final audit starts no new models and creates no full-training-data
model, test predictions or submission.
