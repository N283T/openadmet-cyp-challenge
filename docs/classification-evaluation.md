# Classification candidate evaluation

## Decision to make

Choose both a score recipe and a positive-call rate for each classification
endpoint. The official submission accepts a boolean call for every blind
compound; there is no 49-call requirement. The leaderboard scores only the
subset with a confidently assigned label. A model's broad AUROC/AP, or an MCC
calculated after taking the top 6.53% of *labelled-only* OOF rows, does not
reproduce that decision.

The September 24 [call-count audit](analysis/2026-09-24-d6-call-count-audit.md)
finds different pooled OOF optima across UMAP and scaffold. Treat 49 as the
leaderboard-validated reference count, never as a fixed requirement for a new
model. Select the model and call count jointly outside each evaluated fold;
report instability and do not deploy the pooled OOF maximum as an optimum.
The September 26 [complete-pool rate audit](analysis/2026-09-26-classification-rate-robustness.md)
confirms that the D6 control analog's 49 versus 110 optimum changes by split.
An A4 static-feature proxy also has split-dependent rate selection; it cannot
set the rate for the submitted equal-seven ensemble.

This protocol evaluates candidate generation, call-rate selection and final
deployment separately. It ranks opportunities by estimated MCC improvement and
uncertainty. Diagnostic failures are evidence about the candidate, not an
automatic prohibition on creating one.

## Required prospective prediction bank

For each endpoint and outer split family, save keyed scores for **every** training
molecule, including those with no classification label. The common long schema
is `compound_id`, `endpoint`, `split_id`, `fold`, `model`, `score`, `y_true`;
`y_true` is null for unscored molecules. Every model must cover exactly the same
pool, with one score per molecule from a fit that excluded that molecule's outer
fold. A model that uses auxiliary targets must exclude held-fold molecules across
those targets too. Keep blind-test scores separate.

`evaluate_rates` raises on incomplete pool coverage by default. A
labelled-only sensitivity analysis must pass `allow_partial_pool=True`; its
summary records `deployment_decision_supported=false`. Do not use a diagnostic
partial-pool comparison to select a submitted call vector.

The original saved D6 OOF covers 1,497 labelled molecules, only 24.36% of the
6,146-molecule training pool. The September 24 replay extended the September 15
control and AP-mass90 to the complete UMAP pool and refit both on scaffold folds
([analysis](analysis/2026-09-24-d6-ap-full-pool-replay.md)). This fixes the
population used for retrospective call selection, but neither split reproduces
the observed AP leaderboard loss. The control replay is an analog of the
submitted recipe rather than an exact rerun of historical cross-fitted recipe
selection. Reassigning existing UMAP OOF rows to new folds would leak training
information and is not a new OOF.

The blind set was selected by analog expansion around strong 1A2, 2C9 and 3A4
hits and then measured densely for all four CYPs. A training-only pseudo analog
subset contains 608 molecules but only 164 known D6 TDI labels; its paired
decision interval is too wide to resolve the AP-versus-control leaderboard
direction ([analysis](analysis/2026-09-24-d6-analog-expansion-validation.md)).
Treat this acquisition and label-coverage gap as a limit of the available OOF,
not as something a different target prior or complete-pool call rate has fixed.

## Selection and evaluation

1. Predeclare candidate recipes, split families, call-rate grid and target-prior
   scenarios. Include 49/750 as one historical reference, not a constraint.
2. On each outer fold, choose positive calls from **all molecules in that fold**.
   Compute the score from only those with a supplied classification label; do not
   turn missing labels into negatives. This mirrors the official scoring mask.
3. Select model and call rate using the other outer folds only. Evaluate the
   selected pair on the held-out fold. Compare with the deployed control under
   the same procedure and also show the exact deployed 49-call vector as a
   historical reference.
4. Report the complete rate curve, held-fold choices, paired group-bootstrap MCC
   delta, top-call additions/removals, and results by molecular-support and label
   provenance strata. Report scored and unscored call counts separately.
5. Refit both compared recipes on a prespecified test-like validation family as
   well as UMAP. The
   previous support-matched random/scaffold probes failed on the leaderboard, so
   matching nearest-any Tanimoto alone is insufficient. Inspect label-aware
   neighborhoods and the stability of the choice across families.
6. Apply the chosen recipe and rate to all 750 blind compounds. Review the exact
   changed calls and their uncertainty before submission. The final full-OOF
   optimum is descriptive; only held-fold selection gives an honest estimate of
   the selection procedure.

The estimated hidden D6 positive rate of 7.06% is an approximation derived from
public aggregate leaderboard metrics, not a known test label distribution. Show
sensitivity to plausible priors rather than treating that single value as truth.

## Positive decision evidence

Promote the candidate with the best supported expected gain at the actual
operating boundary. Evidence includes a positive paired held-out MCC estimate,
an acceptable downside under uncertainty, consistent selection of the model and
call rate, and blind additions whose training-label neighborhoods support the
claimed gain. A wide interval or a split disagreement lowers confidence and
can justify a deliberate, limited leaderboard probe when the expected value is
high; it is not by itself a mechanical veto. Record the hypothesis and what the
probe can distinguish before uploading.

Historical leaderboard outcomes can check whether a frozen evaluation protocol
would have separated earlier successes and failures. They are few and selected,
so do not tune the protocol repeatedly against them.

The exact accepted-vector inventory is now recorded in run `4febdf99` and
[summarized here](analysis/2026-09-24-classification-historical-lb-inventory.md).
It fixes seven historical states, the twelve endpoint transitions, and every
changed blind call before an alternative OOF is fitted. Reserve the Sep 18 →
Sep 23 D6 transition as a held-back protocol check.

## Historical leaderboard concordance backtest

Before using a new OOF design to select a future candidate, inventory every
accepted, distinct classification submission with its exact 750-row call vector,
recipe, submission timestamp and official post-scoring endpoint leaderboard
snapshot. Multiple CSV fetches of one leaderboard state are one observation.
The saved 58 classification leaderboard CSVs currently contain seven distinct
N283T `Submitted` timestamps (September 10, 12, 13, 14, 15, 18 and 23), so the
backtest has six chronological transitions, not 58 independent outcomes.
Score D6 and A4 separately; a macro comparison can conceal opposite endpoint
movements. Include the earlier improvements as well as the recent declines.

For each historical recipe, refit the model under the **same new outer-fold
assignments**, rerunning any recipe selection that originally used OOF inside
the training portion of each held fold. Retain the submitted positive-call
count as the first historical replay; later evaluate jointly selected call
rates. Produce scores for the complete held-out molecule pool and evaluate only
the supplied class labels. Never relabel saved UMAP OOF scores as a different
split. Report fold-wise and aggregate estimates with uncertainty.

Compare historical submission *transitions*, not just standalone scores. For
each endpoint transition, record the LB MCC delta and the OOF MCC delta, whether
their signs agree, whether the OOF interval distinguishes the change from zero,
and whether the new OOF ranks the observed winners above the observed losers.
Give special attention to call swaps at the submitted operating points. A useful
OOF design should explain several independent historical directions, including
the September 13 and 15 gains as well as the September 18 and 23 losses. Do not
claim success from matching one dramatic failure.

Freeze the split design and comparison rules before examining its concordance
table; reserve at least one historical transition as a check if the number of
distinct scored submissions permits it. Because all submissions share one
unknown blind set and the public MCC is a bootstrap estimate, even good
historical sign agreement is calibration evidence, not proof of future LB
transfer. Report mismatches and uncertainty rather than adjusting the split
until every sign matches.

## Current implementation and retrospective audit

`classification_decision_evaluation.py` accepts a complete long OOF bank when
available, searches a configured rate grid, selects each model's rate and then
the model/rate pair outside every reported fold, and writes the rate curve,
choices, calls and paired uncertainty. It marks deployment decisions unsupported
when the molecular pool is incomplete.

Recorded run `43665456` applies it to the existing labelled-only D6 control and
AP-mass90 OOF. The AP-mass90 rate selected on other folds varies from 2% to 25%;
the nested candidate-minus-control projected MCC estimate is +0.0211 with a
group-bootstrap 95% interval of [-0.0339, 0.0776]. The all-pool coverage is only
24.36%, and the known leaderboard outcome is negative. This retrospective result
shows why rate search on the labelled-only OOF cannot yet select a submission.
