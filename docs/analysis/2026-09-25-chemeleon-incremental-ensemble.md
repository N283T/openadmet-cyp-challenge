# Adding intermediate CheMeleon readouts to the submitted ensembles

## Question and design

The previous depth screen found complementary predictions among frozen
CheMeleon message-passing steps 1, 3 and 6. This screen asks whether those
readouts improve the **currently submitted multi-family models**. All additions
are fixed, untuned OOF blends at weights 5%, 10% and 20%.

Regression submission components use the parent-retained series split, not the
UMAP split used in the first depth screen. Their folds disagree on 5,043 of
6,384 scored endpoint rows. We therefore refit the three LightGBM readouts on
the parent split in run `regression/20260925T130401763912Z-dcf9b1b8` (70
anchor compounds retained in each training fold; 6,384 nonanchor predictions).
The exact submitted 1A2, 2C9 and 3A4 OOF recipes and the two-component D6
recipe were aligned by compound, endpoint and fold. Run
`regression/20260925T130709627410Z-c8f2824a` evaluated direct pIC50-scale
averages using the official ST-RAE function.

For classification, run `classification/20260925T130752164378Z-d40a0d9c`
reproduced the submitted D6 mean-block cross-fitted OOF exactly (AUROC
`.580764`, AP `.296860`, target-prior MCC `.110109`) and exported its keyed
predictions. Run `classification/20260925T130833379469Z-e0b37526` aligned
these with the submitted A4 seven-member OOF and the original UMAP depth OOF.
Scores were converted to within-fold ranks before mixing, so the weights have
the same meaning across probability and ensemble scores. The comparison uses
only labelled molecules: 1,497 D6 and 3,584 A4.

## Results

Lower regression ST-RAE is better. The four-endpoint macro for the submitted
OOF recipes is `.587127`. Every fixed CheMeleon addition tested worsens it;
the least harmful, 5% final-depth, is `.587304`. The three-depth blend at 5%
is `.587656`. The individual 5% final-depth changes are 1A2 `.670833→.670301`,
2C9 `.454711→.455506`, D6 `.828397→.827620`, and 3A4 `.394567→.395790`.
Thus there is no regression-wide case for adding these simple readouts.

Higher classification AP is better. Relative to the submitted-model OOF after
the common fold-rank normalization:

| Endpoint | Base AP | Fixed addition | AP after addition | AUROC change |
| --- | ---: | --- | ---: | ---: |
| D6 TDI | .295056 | 10% depth 1 | .302755 | .580942 → .579676 |
| A4 TDI | .521602 | 10% depth 6 | .526631 | .805709 → .805485 |

The additions raise labelled OOF AP slightly while AUROC is flat or lower.
Depth 1 against current D6 has Spearman `.544`; final depth against current
A4 has `.722`. These are plausible diversification signals, but this screen
tested several depths and weights and then described the best values. They
are exploratory, not untouched validation results.

## Decision

Do not add a depth readout to the current submission yet. Regression has no
macro benefit under this fixed screen. Classification needs a refit that
predicts the **complete 6,146-molecule held-out pool**, including unlabeled
molecules, followed by held-fold positive-call selection and MCC assessment.
The leaderboard evaluates calls, not AP, and prior D6 labelled-only OOF gains
have reversed on the leaderboard. The next focused experiment is a complete-
pool D6 depth-1 addition and A4 final-depth addition with the present
submission as control, followed by scaffold/analog transfer checks. The
current 750-molecule submitted vectors have not been changed.

Compact keyed OOF and metric files are tracked under
`artifacts/tracked/chemeleon-depth-incremental-20260925/`.
