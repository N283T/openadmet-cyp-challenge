# D6 diversity of the reconstructed native models

The source-style four-member mean improves on its individual members, but it
reaches essentially the same D6 ranking accuracy as the previous global ensemble.
It does not resolve the remaining D6 error.

Recorded diagnostic `20260912T112427694676Z-6eb3beec` reads the completed evaluator
`20260912T112404158683Z-0c4a365a/outputs/all-model-oof.parquet`. Every comparison
uses the same 1,457 canonical non-anchor D6 compounds and folds. Predictions are
fixed; no weights or calibration parameters were fitted. Residual means
prediction minus measured pIC50, and SD uses the population convention.

| Model | Bias | RMSE | MAE | Spearman |
|---|---:|---:|---:|---:|
| Previous global | +0.0282 | 0.8085 | 0.5541 | 0.472752 |
| Four-member mean | +0.0226 | 0.8143 | 0.5508 | 0.472811 |
| Union18 | −0.0078 | 0.8382 | 0.5895 | 0.431863 |
| Isoform6 | +0.0635 | 0.8348 | 0.5712 | 0.415397 |
| Auxiliary8 | +0.0529 | 0.8414 | 0.5804 | 0.411415 |
| Single1 | −0.0182 | 0.8582 | 0.5850 | 0.361932 |

Individual-member prediction Pearson correlations range from **0.560 to 0.707**;
their residual correlations range from **0.876 to 0.942**. The predictions differ
enough to benefit from averaging, while the harder errors remain strongly shared.
The fixed mean therefore improves both RMSE and MAE relative to every constituent.
This is a within-cohort description of averaging, not proof of a gain on unseen
challenge labels.

Against the previous global ensemble, the new mean has prediction correlation
**0.8732** and residual correlation **0.9776**. Its Spearman change is only
**+0.000060**, MAE improves by approximately0.0033, and RMSE worsens by0.0058.
The evaluator's official-reference ST-RAE changes from0.825436 to0.823983;
that tiny endpoint improvement does not establish a ranking breakthrough. Overall
four-target macro ST-RAE worsens from0.582816 to0.622574, so the source mean is
not a replacement for the existing full ensemble.

Biases are small on the measured cohort: the four-member mean's+0.0226 and previous
global's+0.0282 are both much smaller than their approximately0.81 RMSE. A large
common additive offset is therefore not the principal measured-cohort error.
This statement concerns this OOF cohort only.

The blog's reported Spearman0.503 is about0.0302 higher than either ensemble here.
That is a numerical reference across unmatched evaluation cohorts, not a paired
estimate of how much the reconstructed training differs from the author's model.

Outputs include per-model error summaries, prediction/residual correlation
matrices and all fixed equal-pair diagnostics. The latter verify the identity
`MSE((a+b)/2) = (MSE(a)+MSE(b))/2 - mean((a-b)^2)/4`; they quantify disagreement
benefits without searching ensemble weights. The scripts are
`scripts/scp_diversity.py` and `scripts/scp_launch_diversity.py`.
