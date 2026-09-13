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

## Removing nested Caruana entirely

The first pruning round treated nested Caruana as one indivisible component. Its
inner selection is not small: each endpoint/fold retains 4–11 distinct candidate
types and 23–50 greedy additions. Full-data deployment unions the five fold
selections and refits 20/18/19/20 LightGBM models for 1A2/2C9/2D6/3A4, **77 in
total**.

Run `20260912T151747998390Z-f24622ac` removes that component from the candidate
pool and repeats the same target-wise held-out selection over the remaining four.

| Endpoint | Selected count | Cross-fitted ST-RAE | Components |
|---|---:|---:|---|
| CYP1A2 | 3 | 0.687942 | Jeremy frozen TabPFN; Jeremy-score/log2fc TabPFN; Jeremy-finetuned/log2fc LGBM |
| CYP2C9 | 3 | 0.463103 | Jeremy-score/log2fc TabPFN; Jeremy-finetuned/log2fc LGBM; CheMeleon/log2fc TabICL |
| CYP2D6 | 2 | 0.828397 | Jeremy-score/log2fc TabPFN; Jeremy-finetuned/log2fc LGBM |
| CYP3A4 | 2 | 0.411915 | Jeremy-score/log2fc TabPFN; Jeremy-finetuned/log2fc LGBM |
| Macro | 10 endpoint-components | **0.597839** | target-wise equal means |

The no-Caruana macro is 0.001425 higher than the submitted five-component mean
and 0.002240 higher than the first target-wise candidate. The selected fixed
subsets are stable enough that their full-OOF and cross-fitted scores are equal.
This supports removing the 77-model Caruana deployment dependency when management
simplicity is the objective.

Nested selection did have a valid evaluation purpose: the inner group holdout
selected Caruana weights while the untouched outer fold measured that selection.
Without the nesting, choosing weights and reporting them on the same OOF labels
would be optimistic. Nesting is not required at inference, and this comparison
does not find enough incremental held-out value to retain its deployment cost.

Packaging run `20260912T151847989358Z-14af08d1` creates the no-Caruana candidate.
Independent audit `20260912T151916784159Z-8e7cae8c` verifies the selection rule,
all OOF scores and every test mean. Candidate `20260912T151930272078Z-a62d5643`
and prepared submission `20260912T151934226891Z-162bf482` pass dry run. Nothing
was uploaded. This supersedes the earlier 13-component candidate as the preferred
management-oriented pruning result.

## Compressing CYP3A4 Caruana internally

The all-or-nothing removal above did not test the requested middle ground. Run
`20260912T152832952478Z-3e4de7c8` reconstructs every saved CYP3A4 Caruana outer
prediction from its individual LightGBM files (maximum absolute error
1.8e-15), then compares top-K truncation. The deployment-style comparison uses
only the other four folds' Caruana weights to rank candidates for each held-out
fold, refits those candidates on the outer-training labels and evaluates the
unchanged five-component mean.

| CYP3A4 Caruana recipe | Models | Caruana-only ST-RAE | Five-component ST-RAE | Change vs original five |
|---|---:|---:|---:|---:|
| Original nested OOF | up to 11/fold; 20 full-data | 0.456743 | 0.403800 | — |
| Cross-fold top 2, weighted | 2 | 0.556305 | **0.398316** | **-0.005484** |
| Cross-fold top 2, equal | 2 | 0.558623 | 0.398408 | -0.005392 |
| Cross-fold top 7, weighted | 7 | 0.458384 | 0.401037 | -0.002764 |
| Cross-fold top 10, weighted | 10 | 0.452577 | 0.402334 | -0.001466 |
| Cross-fold top 20, weighted | 20 | 0.450564 | 0.403673 | -0.000127 |

The paired conditional bootstrap interval for the weighted top-two change is
[-0.008544, -0.002492]. It improves four of five outer folds; fold 1 worsens by
0.001808. The worse Caruana-only score does not contradict the better final
score: the two retained quantile models complement the other four top-level
components better than the larger Caruana average.

Full-data run `20260912T153211990718Z-e29d7578` replaces CYP3A4 Caruana with a
simple weighted sum of `class_quantile_0.65` (0.520339 normalized weight) and
`morgan_quantile_0.65` (0.479661). This reduces CYP3A4 from 20 LightGBMs to two
and the four-endpoint Caruana total from 77 to 59. Other endpoint predictions
are unchanged.

Packaging run `20260912T153308257696Z-cea40f16` combines this component with the
previous target-wise 3/3/2/5 top-level recipe. Audit
`20260912T153454738585Z-d0cfe511` reconstructs every test mean, verifies exact
equality for the other three endpoints and reports cross-fitted endpoint scores
0.687095/0.463103/0.828397/0.398316, macro **0.594228**. Candidate
`20260912T153532099223Z-96ddb745` and prepared submission
`20260912T153535785660Z-e08fc07c` pass the official dry run. Nothing was
uploaded. This internally compressed candidate supersedes both earlier pruning
candidates, but still retains 20 CYP1A2 Caruana models and is therefore an
intermediate result.

## CYP1A2 internal compression and final operational choice

The retained CYP1A2 Caruana improved the honest top-level component-count
procedure by only 0.000847 versus removing it. Runs
`20260912T162358570573Z-55e7c2fd` and
`20260912T162640128007Z-5d378ca2` therefore repeat the internal top-K comparison
for CYP1A2. The first run evaluates the submitted five-component mean. The
second, decision-relevant run evaluates the deployed fixed three-component mean:
Caruana, Jeremy frozen TabPFN and Jeremy-score/log2fc TabPFN.

| CYP1A2 Caruana recipe | Models | Three-component ST-RAE | Change vs original 20 |
|---|---:|---:|---:|
| Original | 20 | 0.683694 | — |
| Cross-fold top 2, weighted | 2 | 0.692163 | +0.008469 |
| Cross-fold top 4, weighted | 4 | 0.688930 | +0.005236 |
| Cross-fold top 8, weighted | 8 | 0.684808 | +0.001113 |
| Cross-fold top 10, weighted | 10 | 0.683778 | +0.000084 |
| Cross-fold top 17, weighted | 17 | 0.681869 | -0.001825 |

Unlike CYP3A4, CYP1A2 cannot be compressed to two or four models without a clear
point-estimate loss. It needs roughly eight to ten models to recover the original
fixed-ensemble score. That cost is not justified by the 0.000847 advantage in
the honest top-level selection comparison, so the operational choice removes
CYP1A2 Caruana entirely.

Final package `20260912T163036856426Z-bf2e35dc` uses no Caruana for CYP1A2,
CYP2C9 or CYP2D6 and only the two retained models for CYP3A4. Independent audit
`20260912T163050754512Z-ba9c19df` verifies the three non-CYP3A4 predictions are
exactly unchanged from the no-Caruana package, reconstructs all test means and
reports cross-fitted endpoint scores 0.687942/0.463103/0.828397/0.398316, macro
**0.594440**. This is 0.000212 above the intermediate 22-Caruana-model candidate
and uses only two Caruana models rather than 22. Candidate
`20260912T163058242643Z-4f82b69b` and prepared submission
`20260912T163102549701Z-9545d109` pass dry run. The first live attempt was blocked
with 2:30:13 remaining in the confirmed 12-hour cooldown. Retry
`20260913T055138710137Z-e0ce27ff` was accepted at 2026-09-13 05:51 UTC and is
now scored. The MA leaderboard immediately before reflection retained the prior
submission at ST-RAE 0.721167 (rank 67). The reflected candidate scores endpoint
ST-RAE 0.685494/0.369556/1.353707/0.468640 for CYP1A2/CYP2C9/CYP2D6/CYP3A4 and
macro **0.719349** (rank 67), an improvement of 0.001817. This is the preferred
simplified regression candidate.
