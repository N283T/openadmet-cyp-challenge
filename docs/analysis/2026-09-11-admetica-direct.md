# Unchanged Admetica CYP model application

Run `20260911T145352053201Z-cb9e5d6a` applies the four standard Admetica CYP
inhibitor models directly, including their original classification heads.
No challenge fitting, calibration, encoder replacement or TabPFN readout occurs.
Inputs are the current 6,896 standardized molecular structures. Existing
standardization therefore remains in effect; this is unchanged model inference,
not a switch to raw unstandardized strings.

Model assets come from Admetica commit
`64eacdb7aab5a91c83b6f114379455deb3744ee5`, under `ADMET/metabolism/`.
Encoder, predictor and BatchNorm state dictionaries match strictly. CYP3A4 uses
BatchNorm; the other three do not. All modules run in evaluation mode.
Both the unmodified binary probability and pre-sigmoid logit are retained,
and sigmoid(logit) is checked against the predictor's forward output.
All 27,584 probabilities are finite and within [0,1]. All 750 test molecules
have four keyed outputs in `test-probabilities.parquet` (3,000 rows).

## Association with continuous challenge pIC50

For comparison with current models, the table uses exactly the same non-anchor
OOF evaluation molecules. These are direct external pretrained predictions,
not predictions trained by cross-validation in this experiment.

| Endpoint | n | Probability vs pIC50 Pearson | Probability vs pIC50 Spearman | Probability vs FP+Mordred prediction Pearson |
|---|---:|---:|---:|---:|
| CYP1A2 | 1372 | .414507 | .422950 | .500270 |
| CYP2C9 | 1255 | .425625 | .437027 | .415157 |
| CYP2D6 | 1457 | .287561 | .332296 | .490732 |
| CYP3A4 | 2300 | .568048 | .593393 | .574828 |

Logit vs pIC50 Pearson is .415499/.457799/.308865/.594518. Logit ranks are
essentially identical to probability ranks (tiny ties from sigmoid saturation
cause a 2e-6 Spearman difference for 1A2).

On all labeled molecules including anchors, probability/pIC50 Spearman is
.430014/.449820/.338083/.599930. On the neighbor-only subset it is
.339850/.366019/.344620/.663708. All subsets and correlations with the three
existing component predictions are saved in `metrics.csv`.

The unchanged external models are weaker rankers than current challenge-trained
components on matched evaluation rows, but less correlated with FP+Mordred than
vanilla CheMeleon/TabPFN was (.801–.916). This supports retaining them for a
future complementarity test; lower correlation alone does not demonstrate an
ensemble improvement.

The raw outputs are inhibitor probabilities, not pIC50 predictions. Consequently
no pIC50 MAE or ST-RAE is reported for the unconverted output and no regression
submission is prepared. External training overlap has not been audited here;
this run establishes direct application and descriptive associations.
