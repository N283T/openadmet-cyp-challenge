# CYP2D6: what our existing scored submissions identify

Recorded diagnostic `20260912T101837868777Z-b55d06f4` reads only the three existing
upload-attempt prediction files, their saved intents, and matching CYP2D6 leaderboard
CSV snapshots. No network call, new submission, hidden-label query or calibration
constant extraction was performed. Uploaded file SHA-256 hashes match all intents;
all 750 molecule IDs/SMILES align; every N283T leaderboard timestamp matches its
submission exactly.

| Submission (UTC) | Test mean | Test SD | LB ST-RAE | LB Spearman | LB R² |
|---|---:|---:|---:|---:|---:|
| Sep 10 13:32, baseline v1 | 4.592663 | 0.312105 | 1.416195 | 0.373798 | -0.873636 |
| Sep 11 13:27, TabICL replacement | 4.629270 | 0.316658 | 1.412523 | 0.413394 | -0.866707 |
| Sep 12 07:48, full-data ensemble | 4.638156 | 0.405434 | 1.334970 | 0.433905 | -0.792876 |

Means and sample SDs (ddof=1) describe all 750 predictions, **not the unknown scored
half**. The previous full-ensemble report used population SD, explaining its
0.405163 rather than 0.405434; the predictions are identical.

Our D6 ranking improved, while all three models retained similar prediction
centres. The latest model's LB ST-RAE is 0.189774 worse than official TabICL
(1.334970 versus 1.145196), despite Spearman being 0.076750 higher (0.433905 versus
0.357156). This is a concrete discrepancy between ranking and absolute-error
performance. It is compatible with placement problems but does not isolate their
size or direction from these aggregate metrics.

## These are not three placements of one vector

Least-squares affine reconstruction on the same 750 test compounds gives:

| Source → target | Pearson | Affine residual RMSE | Maximum absolute residual |
|---|---:|---:|---:|
| Baseline → TabICL replacement | 0.857821 | 0.162637 | 0.513237 |
| Baseline → full ensemble | 0.857292 | 0.208591 | 0.896894 |
| TabICL replacement → full ensemble | 0.805927 | 0.239862 | 0.985253 |

All pairs fail an exact affine check by substantial margins. Their prediction
ranks also differ. These are model changes, not the fixed-vector affine probes
used by the community code discussed in [the width-calibration report](2026-09-12-d6-interval-width.md).

For a single scored population, each reported R² would supply the moment equation

`R²_i * V = 2 Cov(y,p_i) - Var(p_i) - (mean(y) - mean(p_i))²`.

Even if the scored-half prediction means and variances were known, our three
different models introduce **three distinct unknown covariances**, in addition
to the unknown label mean and variance: three R² equations for at least five
unknown quantities. In the affine-probe construction those covariances are known
multiples of a single covariance, which is the missing simplification here.
Spearman does not supply Pearson covariance. MAE and ST-RAE do not close this
moment system: they depend on the full error/credible-interval distribution.
Furthermore, we do not know which test compounds form the scored half, so the
prediction moments entering these equations are themselves not the saved full-750
moments shown above.

**Conclusion:** the existing three submissions cannot support the blog's blind
moment solve. No blind mean, variance or optimal shift is estimated in this run.
The useful evidence is that model changes improved ordering without substantially
changing the D6 prediction centre; learned offline corrections remain separately
testable, as the nested interval-width experiment demonstrates.

## Evidence and artifacts

`submission-comparison.csv`, `affine-checks.csv`, keyed
`aligned-test-predictions.parquet` and `identifiability.json` are in the diagnostic
outputs. Source snapshots are copied into its declared inputs:

- Attempt `aeaa96e1`, submission `4766cd9f`, D6 snapshot `88b44c6e`.
- Attempt `82910c6b`, submission `f3149927`, D6 snapshot `88b2306a`.
- Attempt `3d2e5c50`, submission `eaafa860`, D6 snapshot `5b722cbb`.

Initial diagnostic `20260912T101754679210Z-6a6ae05a` is retained. Its broad username
filter matched multiple TabICL-containing usernames and therefore omitted the
optional official comparison columns; the final run uses exact `TabICL-baseline`
matching and double-precision descriptive calculations. The receipt/affine
conclusions were unchanged.
