# Diverse CYP2D6 classification models

## Conclusion

A low-capacity Morgan-fingerprint logistic regression adds useful CYP2D6 ranking
diversity to the current ensemble. The preferred candidate uses a 2:1 rank blend
of the current D6 model and logistic regression (`C=0.05`), then calls the top 49
of 750 compounds positive. It replaces 14 of the current 43 calls and adds 20 new
calls. Combined with the 249-call Caruana CYP3A4 candidate, it is a materially new
submission vector rather than an exact endpoint rollback.

The candidate is validated and saved, but was not uploaded. Its improvement is
still an OOF estimate selected from a model screen, so it is not treated as known
leaderboard gain.

## Screen

Run `20260913T112957931554Z-32d00f6e` trains only CYP2D6 on the fixed five UMAP
folds. It compares four genuinely different model families:

- CatBoost on Morgan, Mordred, Morgan+Mordred and frozen external embeddings.
- ExtraTrees on Morgan fingerprints with two leaf-size settings.
- L2 logistic regression on Morgan fingerprints at three regularization levels.
- Squared-similarity-weighted Tanimoto k-nearest neighbors at k=5, 15 and 50.

Every model produces outer-fold OOF predictions and five test predictions. Scores
are rank-normalized within fold before blending with the current model, preventing
fold-specific probability scales from determining the result. Positive counts are
fixed at 43/750 (current count), 49/750 (public top-10 median call rate), or the
inferred 7.06% hidden prevalence.

The more complex models are not best at the decision boundary. Standalone
ExtraTrees reaches AUROC 0.59886 but projected target-prior MCC only 0.03790 at the
43-call rate. CatBoost variants are weaker. The useful signal comes from weak
linear or local-similarity rankers blended with the current ensemble.

## Selected comparisons

Run `20260913T113512499847Z-21e11973` performs a 10,000-resample UMAP-group
conditional bootstrap and packages the candidates.

| D6 candidate | test calls | projected target MCC | delta vs current-43 | 95% conditional interval |
|---|---:|---:|---:|---:|
| Current ranking | 43 | .05717 | — | — |
| 1:1 current/logistic rank blend | 43 | .09309 | +.03592 | [-.00355, .07446] |
| 2:1 current/logistic rank blend | 49 | .09485 | +.03767 | [.00086, .07489] |
| 2:1 current/Tanimoto-k5 blend | 43 | .08777 | +.03060 | [-.00640, .07689] |
| 2:1 current/Tanimoto-k5 blend | 53 | .09668 | +.03951 | [-.00078, .08719] |

The logistic-49 candidate improves projected MCC in all five outer folds relative
to current-43: +.0258, +.0568, +.0066, +.0546 and +.0443. Its interval is only
narrowly above zero and does not correct for choosing it after screening multiple
models and rates. It is therefore promising rather than confirmed.

On hidden test calls, logistic-49 shares 29 positives with the current D6 vector,
drops 14, and adds 20. The Tanimoto-53 candidate is more aggressive, dropping 14
and adding 24, without a stronger uncertainty result. Logistic-49 is preferred.

## Prepared combinations

All packaged files have 750 rows, the official classification schema, and valid
boolean outputs. The primary exploratory candidate is
`logistic_blend_49_caruana_a4.parquet`: 49 D6 positives and 249 A4 positives. Both
endpoints change relative to the current live submission. A lower-risk companion
using the exact old A4 ranking is also retained, along with 43-call logistic and
Tanimoto alternatives.

No upload or submission preparation was performed. If a single classification LB
probe is available and the goal is to avoid repeating an already scored endpoint,
the logistic-49 plus Caruana-A4 candidate is the most informative current choice.

## Submission update

After the analysis was completed, the user explicitly approved this candidate for
submission. Prepared submission `20260913T113807374506Z-77f4825a`, SHA256
`893b79ec83fb19a20676a43cbedfb753ae3388bbe39c38a0293e3f22f2086f6c`, was
accepted at 2026-09-13 11:38 UTC. Automatic pre-upload snapshots saved the MA,
D6 and A4 leaderboards. The first post-upload retrieval still showed the prior
September 12 result. A second retrieval through 11:40 UTC also remained unchanged,
so scoring was not yet reflected at those checks.

The 11:59 UTC retrieval reflects the new result. Macro MCC is **0.352961** at
rank 13, improving from 0.275867 at rank 48. CYP2D6 improves from 0.176698 to
**0.270958** (rank 44 to 21), and CYP3A4 improves from 0.375036 to **0.434963**
(rank 46 to 15). Thus both of the deliberately changed endpoint rankings improve
on the same live submission, and the macro result exceeds the public TabICL
baseline of 0.324967.
