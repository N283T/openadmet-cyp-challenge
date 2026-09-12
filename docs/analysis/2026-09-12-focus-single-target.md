# Single-target challenge adaptation controls

The pure single-target controls did not outperform the submitted ensemble as
standalone models. CYP2D6's adapted representation adds a small complementary
signal: its fixed 25% and 50% blends reach 0.829870 and 0.829681 versus 0.834741.
CYP1A2 does not benefit from either fixed blend. These are OOF results; no new
submission was made.

## Design

Each target independently starts from the same upstream Jeremy medium encoder
bridge `20260912T060426566774Z-7604ffeb`. The challenge adaptation predicts just
one measured pIC50 endpoint. No other target, primary-screen log2FC, Emax, TDI,
or threshold-derived classification label enters encoder adaptation.

The 600-dimensional, depth-4, norm-100 encoder receives a fresh 256-hidden
regression head. The objective is confidence-interval distance normalized by
the inner-training mean-predictor denominator. Encoder/head learning rates are
1e-5/1e-4, batch size 64, seed 42+fold, maximum 160 epochs, patience 20.
Checkpoint selection uses the target's official ST-RAE on a group-separated
10% inner validation set. Label scaling and denominator use inner training only.
The pure-target inner split draws from compounds with that target label, so it
is not necessarily the same inner set as the multitask experiments.

The fixed parent CV `20260911T104703655337Z-a3923255` is unchanged. All labels of
outer nonanchor compounds are excluded from fitting and checkpoint selection;
anchors stay in training. Saved embeddings include all compounds but every fold
uses its own fitted encoder. Independent audits verify the outer membership
exclusions and exact molecule/endpoint/fold/value/confidence-bound alignment.

Readouts use each target's adapted 600-dimensional embeddings plus the existing
four safely nested predicted log2FC features from shared run
`20260912T052549384296Z-bac5fb3a`. The readout is fixed LGBM L1, 300 trees,
15 leaves, minimum 20 samples/leaf, learning rate .05, L2=1. Its training uses
all outer-training target labels. Thus these readouts do use predicted auxiliary
information, while the encoder adaptation itself remains pure single-target.

## Results

| Model | CYP1A2 (1,372 labels) | CYP2D6 (1,457 labels) |
|---|---:|---:|
| Submitted five-model ensemble OOF | 0.683335 | 0.834741 |
| Pure-target direct head | 0.754387 | 0.869362 |
| Pure-target emb600 + nested log2FC → fixed LGBM | 0.718671 | 0.844070 |
| 25% new readout + 75% submitted ensemble | 0.684882 | 0.829870 |
| 50% new readout + 50% submitted ensemble | 0.691269 | 0.829681 |

The two blend weights were prescribed before the readout results. No continuous
OOF weight optimization was performed. The small CYP2D6 blend improvement is a
candidate signal, not evidence that single-task learning solves the weak endpoint.

## Runs and assets

| Target | Encoder adaptation | Direct audit | LGBM readout | Readout/blend audit |
|---|---|---|---|---|
| CYP1A2 | `20260912T082204145061Z-26558f40` | `20260912T082347285554Z-74571659` | `20260912T082600563304Z-a5955327` | `20260912T082616951155Z-ad93e868` |
| CYP2D6 | `20260912T082308319207Z-5c3b657c` | `20260912T082452963642Z-31139312` | `20260912T082618415537Z-79e86a45` | `20260912T082640997403Z-ce95c310` |

Each adaptation saves five fold directories containing best/last checkpoints,
history, memberships, all-compound embeddings and direct predictions. Each
readout saves five fitted LGBM models and keyed OOF predictions. Audits save
`metrics.csv` with whole-target, per-fold and neighbor scores, plus
`verification.json`.

Sources: `scripts/focus_single_target.py`, `scripts/focus_single_readout.py`,
`scripts/audit_focus_single_target.py`, corresponding `launch_focus_single*`
scripts, and immutable configs under `configs/focus-single-target/`.

## Retained failures and scope

- Direct audit `20260912T082330055443Z-64c2fb06` failed because NumPy's `isfinite`
  was called directly on a Polars Series. Conversion to NumPy corrected the
  audit; model fitting was unaffected. The original failure is retained.
- Matching-logFC two-head adaptation `20260912T082451768973Z-fd873af2` was
  interrupted during CYP1A2 fitting following coordination to bound this arm.
  Its record and partial artifacts remain. It supplies no completed comparison
  and must not be interpreted as a negative modeling result. CYP2D6 of that arm
  was never launched. The optional matching-logFC code is retained for future use.
- No additional threshold-derived classification experiment was run in this
  bounded control. Classification is not ruled out by these results.
