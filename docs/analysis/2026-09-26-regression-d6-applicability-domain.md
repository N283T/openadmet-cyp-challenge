# CYP2D6 applicability-domain audit

## Question and sources

Is the D6 regression gap plausibly explained by test molecules being far outside
the labelled training chemistry, and can an applicability-domain (AD) score
identify large D6 prediction errors? The [UNIQUE documentation](https://opensource.nibr.com/UNIQUE/)
defines model-agnostic distance and uncertainty benchmarks. The
[OpenADMET AD article](https://openadmet.ghost.io/applicability-domain-how-to-know-when-your-model-is-actually-relevant/)
emphasizes checking whether a metric actually tracks error; a chemical-distance
threshold alone need not do so. This audit implements corresponding indicators
directly with the project's Polars/RDKit/scikit-learn stack. It does **not** run
the UNIQUE package or reproduce its full error-model pipeline.

Recorded run `regression/20260926T033904883626Z-fd887351` uses 1,493 D6
direct-labelled training compounds and the 750 unique blinded compounds.
For each parent-series or scaffold fold, the D6 training pool contains all
other-fold labels plus retained anchors. Every held-out OOF molecule and every
blind molecule is compared with **that same fold's training pool**, so the
blind comparison does not gain an artificial neighbor-count advantage. The
same 750 blind molecules are evaluated five times under different training
pools; their per-fold values must not be treated as 3,750 independent compounds.

Indicators are nearest-neighbor Morgan radius-2/2048-bit Tanimoto
dissimilarity (top-1 and mean top-5), mean top-5 Euclidean distance in a
training-standardized, training-fitted 32-PC Jeremy medium embedding, and
the absolute z-score of the OOF prediction relative to each fold's D6
training labels. For the parent submitted two-member ensemble, absolute member
disagreement is also available. The error is the per-compound soft distance
outside the official confidence interval. A fold-excluded Random Forest
combines top-5 Tanimoto, latent distance and prediction z-score to predict
that error. It is a diagnostic error-ranking model, not a calibrated prediction
interval or a revised potency model.

The parent error table uses the current submitted two-member D6 OOF control.
The scaffold sensitivity table uses the matched Jeremy-anchor LightGBM OOF;
their error magnitudes are not interchangeable. Association directions can be
checked on both splits.

## Is the blind set structurally outside the labelled D6 training domain?

Lower distance means closer to the fold's D6 training pool. Entries below are
the averages of five fold medians or fractions, not pooled-molecule medians.

| Split | Indicator | OOF median | Blind median | Blind above OOF 90th percentile |
| --- | --- | ---: | ---: | ---: |
| Parent | Top-1 Tanimoto dissimilarity | .675 | **.575** | 4.4% |
| Parent | Top-5 Tanimoto dissimilarity | .713 | **.674** | 4.2% |
| Parent | Latent top-5 distance | 12.81 | **12.46** | 8.9% |
| Scaffold | Top-1 Tanimoto dissimilarity | .659 | **.568** | 5.6% |
| Scaffold | Top-5 Tanimoto dissimilarity | .705 | **.673** | 4.9% |
| Scaffold | Latent top-5 distance | 12.78 | **12.43** | 8.9% |

Under both split geometries, blind molecules are generally **closer** to the
labelled D6 training chemistry than the held-out OOF molecules. This is
consistent with the organizer's analog-expansion sampling around known
non-D6 hits. It does not establish that the blind D6 potency distribution or
model error is in-domain; blind D6 values remain unavailable.

## Do the AD indicators identify D6 error on labelled held-out compounds?

Positive Spearman rho would mean that a larger indicator tends to accompany
larger error. The learned error score is evaluated on a fold that supplied no
error labels to its Random Forest fit.

| Indicator vs. soft error | Parent rho | Scaffold rho |
| --- | ---: | ---: |
| Top-1 Tanimoto dissimilarity | −.028 | −.036 |
| Top-5 Tanimoto dissimilarity | **−.050** | **−.061** |
| Latent top-5 distance | +.027 | −.004 |
| Prediction absolute z-score | +.113 | +.055 |
| Parent member disagreement | +.077 | unavailable |
| Fold-excluded learned error score | +.126 | +.064 |

The farthest top-5 Tanimoto fifth has mean soft error `.296` versus `.384` in
the nearest fifth on parent OOF; scaffold gives `.313` versus `.367`. These
opposite-direction point estimates are small associations with wide practical
ambiguity, not evidence that chemical novelty improves D6 prediction. The
latent distance is nearly uncorrelated with error. The learned combination is
slightly more informative than any single metric on parent OOF, but rho `.126`
is too weak to use as a reliable error warning; it is weaker on scaffold.
Absolute point-error correlations tell the same broad story (see the compact
CSV). No AD rule or molecule exclusion is promoted from this screen.

## Interpretation

The tested **structural** OOD explanation is not supported: blind compounds
are closer to D6-labelled chemistry than the fold-held compounds, and
structural distance does not rank labelled D6 errors. This does not rule out
an **assay-selection** or response-distribution shift. Public D6 direct labels
are strongly selected by D6 primary response, while blind compounds were
selected for other CYPs; the
[selection audit](2026-09-11-label-selection.md) found no public D6 direct
labels among 2,823 compounds with primary log2FC at least −0.5. The
[internal primary-response tests](2026-09-26-regression-d6-internal-primary-tests.md)
found substantial information in measured D6 primary response that was lost
when predicting it from structure. Those observations remain more plausible
diagnostic directions than a nearest-neighbor AD cutoff, although no blind
compound-level D6 outcomes are available to prove the cause.

Compact metrics are in
[`artifacts/tracked/regression-d6-applicability-domain-20260926/`](../../artifacts/tracked/regression-d6-applicability-domain-20260926/README.md).
The keyed OOF and blind AD tables stay in the recorded run. No candidate or
submission changes.
