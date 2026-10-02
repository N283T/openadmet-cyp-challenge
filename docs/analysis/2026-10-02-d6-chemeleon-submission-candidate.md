# Fresh CheMeleon readout: selected D6 submission candidate

**Scored outcome:** [Live evaluation](2026-10-02-d6-chemeleon-live-result.md)
gives D6 MCC .117740060, above the immediately previous .104985814 but well
below historical best .308635. The preparation-stage recommendation below is
now historical; observed-label superiority did not establish Live superiority.

**Submission update:** after explicit approval on October 2, the unchanged
prepared payload was sent and remote receipt confirmed. Scoring is pending;
see [the submission record](2026-10-02-d6-chemeleon-submission-receipt.md).
The unsent statements below describe the preparation and comparison stages.

**Subsequent review:** the user identified that this comparison did not refit
the complete existing eleven-atom ensemble under the new evaluation method.
That [whole-recipe comparison](2026-10-02-d6-legacy-whole-nested.md) is now complete:
400 fits and 20 outer folds. The selected recipe's MCC is .14530/.23066 versus
the refitted historical eleven-atom recipe's .09067/.16759 (UMAP/scaffold).
This supports the recommendation below on observed-label decisions. At a matched
49/750 fraction the historical recipe still wins in both splits; missing weak
outcomes and prospective Live ordering remain unresolved. The exact candidate
remains preserved and unsent.

The continuation produces a concrete new classifier, validates its complete
decision procedure, persists and restores all ten fitted readouts, and prepares
the exact submission payload. The selected candidate is an improvement-seeking
attempt with material transfer uncertainty. It is not a claim that the historical
Live failure has been solved or that sparse development labels predict Live
reliably. No upload has occurred.

## Selected model and decision

The molecular input is a saved 2,048-dimensional generic CheMeleon representation,
not a historical CYP classifier's prediction or ensemble score. Its upstream
checkpoint is a public 47-head ADME pretrain; the original corpus manifest is
unknown. This provenance limits claims about upstream data independence.
Each CYP readout starts with new random parameters.

The readout is 2048 -> 256 ReLU, dropout .15 -> 64 ReLU -> 20 outputs, AdamW
learning rate .002, weight decay .001, 20 full-batch epochs. The D6 class uses
masked binary cross entropy. The four measured CYP primary log2FC targets have
training-only standardized squared-error supervision with weight .2, averaged
over four targets. Other class and paired-dose-response outputs receive zero
loss weight. Unknown TDI classes never become negatives.

Feature means/scales, target means/scales and all thresholds are fitted inside
each training boundary. Two seeds, 20261002 and 20261017, each fit five molecular
group-excluded bag models. Each bag selects a threshold by mean bootstrap MCC
on three deeper group-held-out folds. Final decisions average the ten
probability-minus-trained-threshold margins and call values above zero.
UMAP is the routine deployment grouping; scaffold is a separate robustness
check. No positive-call count, incumbent D6 score, or old positive set is imposed.

Run `20261002T084852768236Z-ef8d81ed` compares matched Morgan, CheMeleon
class-only, and CheMeleon plus primary supervision. All share original outer
groups, deeper fitting boundaries, two seeds and 20 epochs.
Run `20261002T085604402591Z-f4729cb1` validates the actual ten-model averaged
decision and first reconstructs every individual outer call exactly.

| Fully nested ten-model policy | UMAP AP | UMAP MCC | Scaffold AP | Scaffold MCC |
| --- | ---: | ---: | ---: | ---: |
| Fresh Morgan class readout | .29750 | .12572 | .32914 | .16871 |
| Fresh CheMeleon class readout | .30048 | .13178 | .35400 | .21439 |
| **Selected CheMeleon + primary** | **.30345** | **.14530** | **.35805** | **.23066** |

Against the matched Morgan decision, paired MCC delta is +.01958
[-.04450,.08214] UMAP and +.06195 [-.00173,.12174] scaffold. Against the class-only
CheMeleon decision it is +.01351 / +.01627. These are conditional group bootstrap
intervals, not confidence about unseen weak-primary classes or protocol selection.
The decision gain in both splits, not AP alone, selects this new recipe among the
fresh learner options.

## The population problem remains, and the validation is explicit about it

The earlier A4 selection reversal is also tested with this stronger representation
in `20261002T084853832802Z-34215141`. Full-class CheMeleon gives mean broad-known
MCC .32724 UMAP / .39505 scaffold. Hiding weak/missing-primary classes gives
.12850 / .11227. Adding primary supervision does not recover those hidden
outcomes. This rules out dismissing the earlier reversal as only a weak Morgan
representation effect.

Exact partial-label bounds in `20261002T085410573284Z-e4d60d67` distinguish
unknown-outcome uncertainty from resampling known labels. Bounds are verified
against 768 exhaustive small binary assignments and contain every actual A4
masked-control broad-population outcome. They flag ambiguity rather than claim
to recover a reliable historical D6 ordering.

For the frozen complete-pool September 15 -> September 23 AP comparison, the
known true-positive advantage is only two UMAP / seven scaffold calls. Under
equal positive-call counts, three / eight positives in unobserved removed calls,
with unobserved additions negative, suffice to reverse the overall direction.
These are feasible counterexamples on the development prediction population,
not reconstructed blind labels. Exact feasible MCC-difference bounds are
[-.41121,.41681] UMAP and [-.32336,.36702] scaffold, without assuming missing
outcomes or a prevalence. The prior labelled-only successful transitions cannot
calibrate a complete-population selector because their unlabelled predictions
are absent. Missing D6 weak-region labels remain the identifiability limit.

Historical full-pool comparison `20261002T092303063426Z-39c980c3` is deliberately
retained. At the old 49/750 rank fraction and hypothetical prior .07, the new
CheMeleon-primary score gives projected MCC .04910 UMAP / .09308 scaffold,
versus the frozen historical control analog .10356 / .11106. Thus the new model
is **not** established superior at the old operating fraction. The selected new
trained decision has projected MCC .09545 / .15277 in the same hypothetical
prior scenario. These diagnostics do not identify blind accuracy and are not a
grid for choosing a blind positive budget. The historical recipe is an analog,
not a new nested replay of historical recipe selection.

## Final threshold check before preparing the payload

Run `20261002T091033216489Z-692441d6` keeps the network fits identical and compares
observed-label cuts with cuts selected at priors .03/.07/.12 and mean utility over
those priors. The choice rule is recorded before results: choose mean-prior
utility only if its equal-seed outer decision improves at prior .07 and averaged
over all three priors in **both** structural splits. Blind call counts are not
inputs to this rule.

Review `20261002T092041603569Z-24f59b3e` exactly reproduces the original scores
and reconstructs individual calls before averaging margins. Mean-prior utility
reduces ordinary MCC .14530 -> .11870 UMAP / .23066 -> .20127 scaffold. At the
matching .07 scenario it reduces projected MCC .09545 -> .07809 / .15277 ->
.13712; all three scenarios worsen in both splits. The predeclared rule therefore
retains observed-label cuts. This is the selected 175-call payload; no smaller
count was chosen because it looked closer to 49.

## Exact prepared payload and twelve-hour judgment

- Candidate: `20261002T092342108012Z-a9342bae`.
- Prepared submission: **`20261002T092348406577Z-4ea606c1`**.
- SHA-256: `6b23e212457c0f86556dbb200df008a64fb782d46b09659b762608104fb6228d`.
- 750 validated rows. **D6: 175 positives**, 133 additions / seven removals /
  42 retained versus the last sent October 1 payload (140 changed molecules).
- A4: **249 positives**, zero additions/removals, exact last-sent restoration.
- Against the historical best D6 vector: 132 additions / six removals / 43 retained.
- Classification cooldown expired **2026-10-02 11:13:00.944775 JST**.
  A new accepted upload consumes the next twelve-hour opportunity.
- Dry run: valid, `sent=false`.

Fit run `20261002T090029588080Z-2a20da72` persists all ten own readout weights,
scalers, training membership IDs, learned thresholds and blind fold predictions.
Package run `20261002T090500932999Z-5549841c` reloads every weight file, recomputes
all 750 scores per model, verifies agreement within 1e-7, and exactly reproduces
the final decisions. The prepared payload is copied from this immutable output.
The two individual seed policies call 142 / 206 molecules with 64 changed calls;
averaging them is validated, but does not erase fitting uncertainty.

**Recommendation:** use this as the next new-model submission candidate rather
than the earlier 77-call Morgan candidate or a blend of existing D6 scores. It
has a concrete improvement in the validated decision procedure in both splits,
uses broad measured-primary supervision and deliberately abandons the old
positive budget. Spending one twelve-hour opportunity on it is a reasonable
improvement attempt; it is not a promised gain over the historical Live best.
It tests the complete fresh pipeline, not a clean causal isolation of the old
swaps. This judgment keeps the unfavorable historical-rate comparison and the
coverage limitation visible rather than treating development MCC as a guarantee.

No user approval for this specific payload has been received; no upload is
performed. All nine completed new experiment records and manageable own model
weights are retained with input/config/output checksums. The comparison exporter
initially failed on a nullable Boolean mask (`20261002T091445657279Z-1b9b15ba`);
that error is retained, the numeric finite mask is corrected, and the above
comparison uses the successful rerun without changing fitted models. Source
checks and infrastructure tests pass: 73 passed, eight skipped.
