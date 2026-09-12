# Matched observation model implementation

Training run: `20260912T093601580048Z-ec2b7210`.

This is a constrained surrogate observation model. Native signed TDI Emax and TDI-condition pIC50 imply `Emax * sigmoid(ln(10) * h * (pIC50_TDI + log10(C)))`. A learned positive scale and offset per CYP map this response to processed log2fc. It is not an assertion that log2fc equals a raw inhibition fraction. Direct pIC50 is a separate output sharing the learned representation.

Inputs are the verified observation target table `271e3020`, fixed 600-dimensional external embedding `c060d197`, 47 external scores `53eb675c`, current parent assignments `a3923255`, challenge labels `a7063ccd`, and nested log2fc predictions `d272774c` for downstream readout only. The MLP itself does not consume nested log2fc predictions.

Three prespecified arms share 647 → 256 → 128 → 16 architecture, ReLU, dropout 0.1, AdamW (learning rate 0.001, weight decay 0.0001), minibatches at most256 and seed42+fold:

- `independent`: all16 outputs free (4 direct,4 TDI pIC50,4 Emax,4 log2fc).
- `hill_fixed`: log2fc outputs replaced by the Hill/affine surrogate, h=1.
- `hill_learned`: same relation with one bounded h per CYP, `0.25 + 3.75 sigmoid(hraw)`, initialized at1 and penalized by0.01 mean(log(h)^2).

Training loss averages head-specific observed labels: direct interval distance normalized by the training constant-predictor interval error, weight2; standardized TDI MSE weight0.5; standardized Emax MSE weight0.25; standardized log2fc MSE weight0.5. No absent task labels are imputed into the loss. Direct denominators and target standard deviations have a0.05 numerical floor. The learned positive affine scale uses softplus.

All outer validation compounds are excluded from every target head. The inner split holds out10% of non-anchor groups; every group containing an anchor remains entirely in inner training. Inner and outer phases fit their own feature/target scaling on their respective training membership. Inner direct ST-RAE selects epochs, at most150 with patience20. A fresh network with the original seed then refits all outer training molecules for exactly the selected epoch count; no outer score selects training duration.

Each final checkpoint contains model state, feature column order, x_mean/x_std,16 target means/stds, direct normalization denominators and selected epochs. Full per-epoch inner diagnostics and molecule memberships are saved. Full6896-row128-dimensional latent representations are saved for each outer fold/arm. The fixed LGBM readout uses these128 features plus the same four nested log2fc features and the same300-tree/15-leaf configuration as other comparisons. Both direct and LGBM predictions carry molecule, endpoint, fold and split identifiers.

The serialized16-output layer is kept structurally identical in all arms; its four free log2fc coordinates are unused in curved arms. The affine parameters receive gradients only in the curved arms. This preserves the shared representation and direct/TDI/Emax head architectures across the control and constrained treatments.

## Completed outcomes

| Prediction | Endpoint | Independent | Fixed h=1 | Learned global h |
|---|---|---:|---:|---:|
| Direct MLP | 1A2 | 0.734246 | 0.743712 | 0.737859 |
| Direct MLP | 2C9 | 0.509909 | 0.522679 | 0.515186 |
| Direct MLP | 2D6 | 0.886731 | 0.904206 | 0.894060 |
| Direct MLP | 3A4 | 0.439728 | 0.442683 | 0.441168 |
| Latent + nested LGBM | 1A2 | 0.729005 | 0.749733 | 0.748530 |
| Latent + nested LGBM | 2C9 | 0.544555 | 0.539907 | 0.528820 |
| Latent + nested LGBM | 2D6 | 0.892313 | 0.911272 | 0.902257 |
| Latent + nested LGBM | 3A4 | 0.466992 | 0.470412 | 0.462711 |

ST-RAE, lower is better. The constrained arms do not improve either focus endpoint relative to the matched independent arm. This small-MLP control itself is weaker than the existing strong components; no model is promoted. This is evidence about this fixed implementation, not evidence against a dose-response relationship in the underlying assay.

Selected epoch counts by outer fold0–4: independent `[2,3,3,2,1]`, fixed Hill `[5,3,4,2,1]`, learned Hill `[2,3,4,2,1]`. The fitted global Hill exponents stay near initialization (roughly0.98–1.05). Epoch selection was not overridden after seeing these results.

The run completed all15 neural refits and60 LGBM readouts. `observation-oof.parquet` contains6 methods×6,384 rows; `metrics.csv` contains all24 endpoint summaries; `memberships.parquet` contains five outer-fold assignments for every compound. Detailed training selection scores are in `training-history.parquet` (the summary's `inner_ST_RAE` field is null because refit does not score an inner validation set). Root-level independent audit and joint candidate integration are separate follow-ups.
