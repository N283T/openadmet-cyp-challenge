# Robust one-sided D6 primary-bound follow-up

Compact outputs from the recorded runs in the
[analysis](../../../docs/analysis/2026-09-26-regression-d6-robust-primary-bound-followup.md).
Full keyed OOF, assignments, bootstrap draws, input snapshots and provenance
remain in immutable local runs.

| File | SHA-256 |
| --- | --- |
| `standalone-metrics.csv` | `09b7499266a64f24b065d213470e4b0601544d19b0b8c25089184f59ad88c426` |
| `blend-metrics.csv` | `3f00b6803b7a8cecf1f3a5b4b832932d006179c444ead05afc868d9d85592f82` |
| `fold-metrics.csv` | `2537e01522a157bf81baea1eb0d0b4f8eb537146cadb2d728613dda99827dc8c` |
| `bootstrap-summary.csv` | `c51478dbbb87ce5f1a8996d1bce878c8bd5f91d8226e55f2d8e560d76cce94d3` |
| `synthetic-metrics.csv` | `e8b830eba86ffe9538f9c93e78d87601199343ae75119f20b658580caffb8bac` |
| `lofo-selection.csv` | `2c2d9d835e19d2e25e75efd51147be8a6d689f14a7cd04bf986b732b4f78d23c` |
| `lofo-metrics.csv` | `af06bba0da7c215b401cefb7935d6846bf85eea41c0aa8c677b0dbb656923571` |
| `candidate-comparison.csv` | `018cadef31ca9069e635001162be6f21394666310a7857fc10896631e60a6c08` |
| `fullfit-summary.json` | `19bd9bf04b4af3023254f6bc7f07437626216346a2ec23249fec30c3ba450677` |
| `d6-robust-bound-model.txt` | `a036e12b69569328f02867d7b02e5505d902aeaf7dbd7e95c4ef99fce31d8578` |
| `prepared-predictions.csv` | `399cd07c66c0d4e7e1fb468b6a36a471d1572951e7f10e324e0b597d216c7b29` |

The LightGBM booster predicts a residual centered on the D6 direct-training
mean `4.784202026339585`. Add that mean to its raw output. The feature order
is sorted `score_` columns followed by sorted `emb_` columns. Full deployment
is defined in `scripts/regression_d6_primary_robust_fullfit.py`: match the
new raw D6 vector's mean and SD to the September 24 submitted D6 vector,
then blend equally. Prepared submission `20260926T090548323924Z-46188564`
was later uploaded. Its reflected
[Live D6 result](../../../docs/analysis/2026-09-26-regression-d6-robust-live-result.md)
is ST-RAE `.549341` versus the previous `.565520`.
