# CYP2D6 complete-pool OOF transfer check

## Question and runs

Does including every held-fold molecule in positive-call selection change the
ranking of D6 readouts compared with the historical labelled-only OOF? Runs
`c8368f50` (balanced LightGBM) and `2542e03a` (ExtraTrees) refit two fixed
external-feature readouts on the same UMAP and scaffold assignments. They save
one fold-excluded score per model and each of the 6,146 training molecules;
1,497 have D6 labels and 4,649 retain null labels. Both representations are
frozen, and readout fitting and ExtraTrees feature selection exclude the held
fold. Their labels-only score rows reproduce the separately saved OOF scores to
floating-point precision.

Runs `d67281d0` and `cc83a23c` evaluate the complete pool and the old
labelled-only boundary, respectively comparing two LightGBM feature sets and
LightGBM with an ExtraTrees readout. They prespecify the historical 49/750
rate, a 2–20% rate grid, and scored-population priors of 5%, 7.06% and 10%.
The evaluation selects model and rate without the reported fold and uses paired
group bootstrap intervals. Final compact outputs and complete-pool scores are
promoted to `artifacts/tracked/d6-full-pool-transfer-20260924/`.

## Fixed historical rate

At the working 7.06% prior, the candidate-minus-control projected MCC is:

| Comparison | Split | All-pool boundary | Labelled-only boundary |
|---|---|---:|---:|
| LightGBM with Jeremy scores versus LightGBM without them | UMAP | −.00369 | .00000 |
| Same | Scaffold | +.01470 | −.00521 |
| Jeremy-score ExtraTrees versus static LightGBM | UMAP | +.00327 | −.00498 |
| Same | Scaffold | −.01404 | −.02550 |

Thus the population correction can change the direction of a comparison, as
seen with the scaffold LightGBM pair. It also changes which labelled molecules
receive the fixed call budget: at the full-pool boundary roughly 84–89 of 400
fold-level calls are labelled, versus 97–99 of 97–99 when the unlabelled pool
is excluded. Every fixed-rate paired bootstrap interval includes zero.

## Joint model/rate selection

The selected operating rate varies materially across folds. For the UMAP
complete-pool LightGBM comparison, the control chooses 6–20% and the
Jeremy-score variant 2–6.53%; the pooled candidate-minus-control estimate is
−.0222 with 95% interval [−.0618, +.0163]. On scaffold the corresponding
estimate is −.0006 [−.0321, +.0312].

For the distinct ExtraTrees challenger, UMAP complete-pool selection estimates
−.0065 [−.0445, +.0304], while scaffold estimates +.0465
[+.0067, +.0871]. That positive scaffold result is one of several comparisons
and reverses on UMAP. It is evidence of split sensitivity, not a supported
deployment gain.

## Decision and remaining transfer test

Full-molecule call selection is now executable and auditable for two genuinely
refit model families. It removes a known mismatch between the earlier 1,497-row
OOF and the 750-row blind ranking. It does not by itself make local gains track
the leaderboard: these are bounded static-feature readouts, not exact replicas
of the September 15 control or September 23 AP-mass90 package. The discordant
UMAP and scaffold assessments argue against promoting either tested candidate.

The next necessary experiment is to reproduce the submitted September 15
control and AP-mass90 recipes on the same complete-molecule outer folds, with
fold-local feature/model selection and auxiliary-label exclusions. Evaluate
their historical 49-call boundaries before optimizing rates. Compare the
direction of the September 14 → 15 gain and September 15 → 18 loss, then use
September 18 → 23 as a procedural held-back check. Saved labelled-only scores
cannot be reassigned to new folds without leakage.
