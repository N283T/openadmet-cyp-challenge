# D6 failures, matched refits and a concrete prospective candidate

Subsequent decision: the user requested additional local validation before
consuming a twelve-hour submission opportunity. The
[complete-policy validation](2026-10-02-d6-submission-value-validation.md) now
evaluates the actual five-model algorithm and controlled alternatives. The
prepared candidate below remains preserved and unsent; upload is not recommended
on that additional evidence. Its earlier candidate rationale is retained here.

## Purpose and decision

The user resumed the investigation and clarified that new models and new splits
are welcome, but experiments must connect to earlier failures and proceed to a
concrete decision. The prior disposition closing the investigation is historical.
Creating a network alone did not answer the diagnosis. This work compares actual
historical recipes on matching splits, audits prior failed decisions, and prepares
an independent molecular-model submission candidate with a stated prospective
question. Uncertainty is evidence for that decision, not an automatic veto.

The question is whether an independent D6 molecular learner, train-only decision
thresholds and fold inference transfer better than the failed AP-selected rank
mixture. This submission tests the combined pipeline. It cannot separately identify
the effects of architecture, thresholds, aggregation and changed positive count.

## What connects to the failures

1. October 1 changed identities at a fixed count: 21 D6 calls were added and 21
   removed, with 49 positives in both payloads. Live MCC fell .308635 to .104986,
   and both precision and recall fell. Count alone cannot explain that reversal.
   The retained conditional multi-block AP strategy increased AP .296860 to
   .318808, while ordinary historical-rate decision MCC fell .158187 to .151597
   and the prior-.07 projection fell .109717 to .104588. AP did not establish that
   the substituted positive calls were better. These conditional strategy metrics
   are not an independent validation of the exact full-development-selected
   blind mixture, and do not prove the molecular cause of its Live loss.
2. The earlier AP-mass90 failure is a counterexample to a universal explanation
   based on AP alone, incomplete prediction pools, or one split. September 24
   complete-pool control/AP refits on UMAP and scaffold still favored AP at the
   historical rate, despite its Live loss (.30864 to .23483). Density weighting
   further favored that failed candidate. Those checks already existed; merely
   repeating Morgan or changing splits is not itself a new diagnosis. See the
   [complete-pool replay](2026-09-24-d6-ap-full-pool-replay.md).
3. Prior support6 calls lost label-aware neighborhood support, but AP-mass90 did
   not share that pattern. A universal neighbor gate would not explain both.
   The selected blind analog population and unknown D6 scored-label mask remain
   plausible transfer problems, not established causes. See the
   [three-submission retrospective](2026-09-24-d6-three-submission-retrospective.md).
4. The new matched fits expose a second practical sensitivity: nearly unchanged
   ranking can produce very different positive calls after refitting and threshold
   transfer. For the fresh UMAP network, full-refit versus centered-fold inference
   has rank correlation about .975, but gives 158 versus 77 positives and changes
   81 calls. This is a newly measured mechanism of decision instability; it does
   not establish that it caused any past Live failure.

## Matched models and splits

Training run `classification/20261002T021327768010Z-c37abd9c` refits historical
Morgan logistic C=.05, RBF SVC C=1, ComplementNB alpha=1 and the fresh D6 network.
Morgan models were already tested in September; this is a matched decision audit,
not a claim that Morgan is a new model. Exact original five-fold UMAP/scaffold
assignments are inputs. Three inner grouped folds exclude all measurements of
held-out molecules. Inner score quantiles select bootstrap-mean MCC thresholds;
there is no 49-call constraint. All 6,146 development molecules receive molecular
scores; the 1,497 labelled D6 molecules supply classification supervision.

| Model | UMAP AP | UMAP nested MCC | Scaffold AP | Scaffold nested MCC |
| --- | ---: | ---: | ---: | ---: |
| Historical Morgan logistic | .2710 | .0843 | .3073 | .1086 |
| Historical Morgan RBF | .2758 | .0777 | .3284 | .1246 |
| Historical ComplementNB | .2898 | .0895 | .3303 | .0833 |
| Fresh D6 network | .2939 | .1088 | .3188 | .1698 |

These MCC values are bootstrap means conditional on the retained outer fits.
They do not include training/selection uncertainty. Scaffold ComplementNB has
higher AP than the network and markedly lower decision MCC, illustrating the
objective mismatch without invoking Live labels.

The historical whole ensemble is stronger than these historical atomic Morgan
recipes. Retrospective full-pool rank-rate projections at assumed prior .07 have
grid maxima .1089 versus .0855 on UMAP and .1344 versus .1094 on scaffold for
historical raw ensemble versus fresh network. These maxima are post hoc diagnostics,
not honest performance estimates or known blind priors. The new network does
**not** establish superiority to the best whole ensemble by that proxy. Missing
labels remain missing in all full-pool calculations.

## Why this is a submission candidate

The candidate deliberately tests an independent molecular learner and an explicit
probability decision policy, rather than mixing the old scored vector into another
AP improvement. Routine UMAP is the deployment split; scaffold is a sensitivity
check, not selected because its MCC is larger. For each of five UMAP outer models,
subtract its own inner-selected threshold from its blind score. Average the five
margins equally and call positive when the mean is greater than zero. There are
no fitted mixture weights, old-score inputs, fixed positive count, or hidden-label
reconstruction. The outcome is 77 D6 positives.

This policy was selected after examining the fits and inference sensitivity.
The outer OOF audit evaluates individual held-out fold decisions, not an identical
five-model bag. The final bag therefore lacks a directly matching nested-CV
estimate. It is an explicit prospective probe, not a validated winner. Preparing
it makes the uncertain hypothesis reviewable and testable; the negative whole-
ensemble proxy comparison is retained rather than omitted or made a veto.

If Live improves, this supports transfer of the combined independent pipeline,
not one isolated causal explanation or a proven optimum of 77 calls. If it
worsens, that is an informative failure of this combined probability/fold policy;
it cannot be blamed solely on AP or 49 calls. Preserve the exact model, per-fold
scores, thresholds, alternate inference vectors and payload for either outcome.
An unsent local candidate cannot resolve the unknown blind-label transfer question.

## Prepared payload and validation

Final review run: `classification/20261002T022712675753Z-4f8c8563`.
Preliminary review `20261002T022503807478Z-ec177efd` is preserved and superseded:
it initially audited only the single-model AP strategy; the final report adds
the multi-block strategy linked to the October 1 probe. Payload calls are unchanged.

Candidate: `20261002T023124733954Z-c7da0ec3`.
Prepared submission: `20261002T023446563579Z-e264a4bf`.
Payload SHA-256: `0ad9d281b58403881500405832a9d5047e0e43f02dfd11659d0ae4729e9a89eb`.

| Target | Positives | Added vs last sent | Removed vs last sent |
| --- | ---: | ---: | ---: |
| D6 | 77 | 54 | 26 |
| A4 | 249 | 0 | 0 |

D6 retains 23 calls from the last sent failed probe. Against the best-scored D6
vector it adds 46, removes 18 and retains 31. A4 is the exact restored best-scored
249-call vector. Regression is unchanged. Preparation and dry-run validate all
750 rows and the four required columns. Dry-run returns `sent=false` and next
eligibility October 2, 11:13:00.944775 JST, already elapsed at preparation.
No challenge upload has occurred. The remaining upload requires the user's
specific approval for this exact prepared payload under AGENTS.md.

Source, configs, manageable weights, keyed predictions, diagnostic curves and
negative results are retained under `artifacts/tracked/cyp-d6-zero-start-*` for
the three run IDs above. The final run's generated report and configuration
provide the complete aliased input provenance and comparison bank.
