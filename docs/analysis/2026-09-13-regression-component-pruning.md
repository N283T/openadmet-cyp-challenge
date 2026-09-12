# Target-wise pruning of the submitted regression ensemble

The submitted regression model averages five top-level components for every
endpoint. This analysis tests whether each endpoint needs all five. It reuses the
saved outer-fold predictions and trains no new base model.

## Selection protocol

For each endpoint and requested component count from one through five, enumerate
every equal-weight subset. For each outer fold, select the best subset using only
the other four folds and evaluate it on the held-out fold. Concatenate the five
held-out predictions to score the component-count procedure. The deployment
subset at the chosen count is fitted on all OOF predictions and is reported
separately from the cross-fitted estimate.

The selected count is the smallest whose paired 1,000-replicate conditional
bootstrap interval versus the best cross-fitted count includes zero. This
uncertainty is conditional on saved base predictions and does not include base
model refitting or prior model-family screening.

Selection run: `20260912T150626214420Z-d602ffed`. An earlier completed run
`79f26288` contains the target-wise analysis before the shared-family table was
added; no reported value changed.

## Result

| Endpoint | Five-component ST-RAE | Selected count | Cross-fitted ST-RAE | Fixed deployment subset |
|---|---:|---:|---:|---|
| CYP1A2 | 0.683335 | 3 | 0.687095 | nested Caruana; Jeremy frozen TabPFN; Jeremy-score/log2fc TabPFN |
| CYP2C9 | 0.463781 | 3 | 0.463103 | Jeremy-score/log2fc TabPFN; Jeremy-finetuned/log2fc LGBM; CheMeleon/log2fc TabICL |
| CYP2D6 | 0.834741 | 2 | 0.828397 | Jeremy-score/log2fc TabPFN; Jeremy-finetuned/log2fc LGBM |
| CYP3A4 | 0.403800 | 5 | 0.403800 | all five components |
| Macro | **0.596414** | 13 endpoint-components | **0.595599** | target-wise equal means |

The corresponding full-OOF-fitted fixed subsets score macro **0.594749**. That
number is a deployment-fit diagnostic, not another held-out estimate. The
cross-fitted result is the relevant comparison.

The target-wise recipe reduces 20 endpoint-component uses to 13, or 35%. It does
not eliminate a complete top-level family because CYP3A4 retains all five. A
shared three-family subset—nested Caruana, Jeremy-score/log2fc TabPFN and
Jeremy-finetuned/log2fc LGBM—is selected identically in all five held-out-fold
fits and scores cross-fitted macro 0.597496. It would remove the frozen-medium
TabPFN and CheMeleon/TabICL dependencies at a point difference of 0.001082 from
the five-component mean. This shared-family option is preserved as evidence but
is not the packaged target-wise candidate.

## Packaged candidate and verification

Packaging run `20260912T151057434347Z-ba71e814` averages the saved full-data
component predictions according to the endpoint table. It validates all 6,525
training memberships, zero train/test overlap, all 3,000 long predictions and
the official 750-row schema. The first packaging attempt `8d286e61` failed before
reading model outputs because snapshot inputs that already had numeric prefixes
received a second prefix. The corrected run uses the actual immutable snapshot
names; the failed attempt remains recorded.

Independent audit `20260912T151306121315Z-97a2dcab` verifies the component-count
rule, reconstructs every selected OOF score and every packaged test mean, and
reproduces cross-fitted macro 0.595599 and full-OOF macro 0.594749. Candidate
`20260912T151111373790Z-b45b0614` and prepared submission
`20260912T151118536196Z-ddd9e335` pass the official dry run. Nothing was uploaded.
