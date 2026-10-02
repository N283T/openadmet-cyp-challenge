# D6 recommendation audit: an unresolved validator was mistaken for a model selector

The user's standing request was to connect validation to previous leaderboard
failures and establish what supports an improvement decision. Building a fresh
model and proving a same-method CV difference did not answer that question.
The stated confidence in spending a submission opportunity exceeded the evidence.
The user had allowed new models and splits as investigative tools, not as a
replacement for this diagnosis. This is a failure of the recommendation's
inference, not simply an unlucky model outcome.

## Evidence already available before the recommendation

The September 23 AP-mass90 change lost .07381 Live MCC versus the historical
best. Complete-pool recipe analogs still predicted a gain on both splits.
Run `20261002T041315546883Z-49b9e5f8` had already measured:

| Diagnostic | UMAP estimated MCC change | Scaffold estimated MCC change | Actual Live change |
| --- | ---: | ---: | ---: |
| Uniform, hypothetical prior .07 | +.00708 | +.04536 | -.07381 |
| Primary-selection weighting | +.00078 | +.04288 | -.07381 |

Weighting did not resolve the reversal. Earlier successful changes had matching
signs, but their diagnostic call pools were labelled-only; counting them together
with complete-pool failures as independent validation successes would be invalid.
These facts were explicit in the earlier reports, not discovered after the new
submission. The primary label-coverage gap and the A4 masking reversal likewise
predated the recommendation. Merely retaining caveats alongside a recommendation
did not establish a justified transition from model CV to Live confidence.

## The whole-model comparison also changed the reference decision policy

Audit run `20261002T130706977633Z-daff4e97` uses saved, identical two-seed
outer predictions and verifies the prior rank-fraction diagnostic exactly.
It compares the selected model's trained decisions with two legacy references:
the new inner-selected threshold policy, and the old 49/750 rank-fraction
diagnostic. The latter is a historical reference, not a count constraint on a
new model or a proposed restored payload.

| Split | Legacy new-threshold MCC | Legacy historical-fraction MCC | Selected trained MCC |
| --- | ---: | ---: | ---: |
| UMAP | .09067 | .14674 | .14530 |
| Scaffold | .16759 | .14501 | .23066 |

Against legacy new thresholds, the selected UMAP advantage is +.05462,
conditional 95% interval [.00306, .10955]. Against the historical fraction it
is **-.00145**, interval [-.06786, .06298]: the routine UMAP advantage disappears.
Scaffold still favors the selected model by +.08565 against that fraction,
interval [.01309, .14764]. This is not a universal refutation of the fresh model,
but it invalidates interpreting the earlier two-split result as superiority
over the actual historical model-and-operating-policy combination.

The historical-fraction reference uses freshly nested recipe analogs and ten
readout averages, not an exact refit of every historical discovery/deployment
step. All intervals condition on saved fits and omit selection uncertainty.
This is a post-result comparator audit, not an independently validated selector.
It isolates a decision-policy confound; it does not explain the exact molecular
cause of the Live failure, validate 49 calls, or repair missing class coverage.

## What the completed audit establishes

Two problems in the recommendation are now explicit: unresolved historical
sign reversals were bypassed, and a fair model comparison under new thresholds
was extended into a stronger claim about a previously deployed policy. The
fresh candidate's Live MCC .11774 remains below the historical .30864.
Adding another network, or treating a conditional interval excluding zero as
prospective confidence, does not repair either inferential problem.

The investigative objective remains validation transfer: compare submitted
model-and-decision policies, preserve contradictory historical results, and test
proposed explanations against both failures and successes with explicit source
coverage. New models may be probes for those explanations; their CV ranking
alone is not completion of that objective. Missing D6 outcomes remain missing,
and the A4 masking control supplies a mechanism analogy rather than proof of
the realized D6 error. No replacement selector has been validated by this audit.

No new model is trained, no candidate is created, and no upload is performed.
The audit retains all source-backed contradictions and its complete derived
tables with provenance. The arithmetic/identity assertions and `pixi run check`
pass. This report corrects the recommendation's scope without reinterpreting
the old 49-call payload as the desired modeling baseline.
