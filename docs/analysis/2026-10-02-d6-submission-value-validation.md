# D6: validate the complete decision pipeline before spending a submission

## Decision

The user requested further local validation because a classification submission
consumes a twelve-hour opportunity. The original 77-positive candidate is preserved,
but is not recommended for upload on this evidence. No challenge upload, new blind
payload, incumbent restoration or A4/regression change occurred. This is a decision
based on completed additional experiments, not a reversal in response to the user's
tone or an automatic rejection of uncertainty.

The original candidate's exact five-model decision algorithm now has outer-fold
evaluation. It has not demonstrated a consistent benefit over a single refit,
and its positive decisions remain sensitive to repeated fitting. Two controlled
ablations reduce some numerical problems without establishing a reproducible
replacement. Creating a valid file was not sufficient evidence for a submission.

## Frozen protocol and boundaries

Run `classification/20261002T031830189787Z-23751617` evaluates the original
100-epoch network and three actual historical Morgan recipes. Five original outer
folds are fixed for UMAP and scaffold. Two seeds repeat the inner fitting and bag
assignments. Inside each outer training pool, build a five-model bag; each model
uses three additional grouped folds to select its threshold. All outer test
molecules and their measurements are excluded from every fit and threshold choice.

For every model and seed, retain four prespecified policies:

- **Centered fold mean:** average five score-minus-own-threshold margins, then
  call positive if the average exceeds zero. This is the original candidate policy.
- **Majority vote:** at least three of five own-threshold decisions are positive.
- **Pooled OOF cut:** select a common cut from training-only bag OOF scores and
  apply it to the mean of five predictions.
- **Single refit:** inner-three-fold cut, followed by one model fit on all outer
  training molecules; the comparator for the original earlier model audit.

Every held-out pool includes all development molecules, including those with no
D6 class label. Only supplied labels enter losses, threshold MCC or scored MCC;
missing labels never become negatives. There are 6,146 development molecules and
1,497 known D6 labels. No 49-call budget or assumed blind prevalence enters the
threshold selection. No cache of historical scores enters a new model fit.

This evaluates the complete relative algorithm on a smaller training sample:
each bag model sees 64% of development versus 80% in blind deployment. The single
refit sees 80%. Cached historical whole-ensemble scores use their original outer
fits and are only matched rank-rate diagnostics, not new nested threshold learners.
The outer labels were inspected in earlier experiments; this is conditional
protocol validation, not an independent model-family discovery dataset. The two
seeds change both fitting randomness and nested partition assignments; their
comparison is whole-pipeline reproducibility, not an isolated initialization test.

## Original five-model policy results

The numbers below are ordinary held-out MCC, not the earlier conditional bootstrap
mean. Each cell lists seed 20261002 / seed 20261017.

| Policy, 100 epochs | UMAP MCC | Scaffold MCC |
| --- | ---: | ---: |
| Centered fold mean | .1005 / .0937 | .1765 / .1522 |
| Majority vote | .0999 / .0779 | .1544 / .1498 |
| Pooled OOF cut | .0463 / .0292 | .1217 / .1519 |
| Single refit | .1086 / .1046 | .1708 / .1313 |

Original-policy AP is .2998/.2978 on UMAP and .3212/.3139 on scaffold. There is
ranking signal, but the averaged decision's extra benefit is not established.
Seed-averaged centered-minus-single MCC is -.0095 on UMAP, conditional molecular-
group bootstrap 95% interval [-.0397,.0240]; scaffold +.0133 [-.0150,.0448]. All
comparisons against the three historical atomic Morgan policies also cross zero.
Intervals keep repeated seeds paired and omit refitting and protocol-selection
uncertainty. A zero-crossing interval alone is not the reason to reject an upload.

The original policy calls 15.1% / 20.5% of full UMAP held-out pools and 22.0% /
37.5% of scaffold pools. Across seeds, 367 UMAP and 1,132 scaffold molecule calls
change among 6,146. Positive-set Jaccard is .713/.528. UMAP score correlation is
.9948 despite the decision swaps. These held-out percentages are not forecasts
of the 750-blind call rate or constraints on it.

## Controlled duration ablation

Run `classification/20261002T032805271125Z-bdf341bd` changes only network training
duration from 100 to 20 epochs, retaining architecture, splits, seeds, inner
selection and the four policies. This is one prespecified follow-up after the
first run, not a duration grid searched for a winning blind payload.

| Centered fold mean | UMAP MCC | Scaffold MCC |
| --- | ---: | ---: |
| 100 epochs | .1005 / .0937 | .1765 / .1522 |
| 20 epochs | .1377 / .1066 | .1649 / .1450 |

Shorter fitting reduces the prediction-scale shift. Median refit-validation
90th-percentile score minus inner-validation 90th-percentile score falls from
.2473 to .0264 on UMAP and .1628 to .0303 on scaffold. The populations differ,
so this is a scale diagnostic, not a same-molecule causal decomposition.
Threshold ranges contract from approximately zero-to-one at 100 epochs to
.0023-.4565 UMAP / .0114-.4703 scaffold at 20 epochs.

Mean-prediction Brier errors fall from about .204/.203 to .174/.169 on UMAP/
scaffold. A pooled empirical constant-rate reference is .1696; it is an
illustrative descriptive baseline, not an independently trained competitor.
AP remains about .296-.298 UMAP and .326-.330 scaffold. Less extreme scores
therefore do not automatically give a stronger ranking or decision rule.

Mean MCC change for 20-minus-100 is +.0250 UMAP, conditional interval
[-.0267,.0701], and -.0094 scaffold [-.0488,.0286]. Call identity reproducibility
does not improve generally: centered-policy changes across seeds are 1,040 / 1,004
and positive Jaccard .635/.614. The corrected score scale does not by itself solve
threshold selection and deployment.

## Controlled rank-transfer ablation

Run `classification/20261002T033604148307Z-17b6eddf` returns to the original
100-epoch fits and adds one inference policy. For each model, convert outer test
scores to percentiles within its complete held-out pool, and convert the selected
inner cut to its percentile among all inner-OOF outer-training molecules. Average
the five percentile margins and call above zero. No test labels enter ranking or
the selected percentile, and no historical call count is prescribed.

All four original policy score and call vectors reproduce the first run exactly;
the comparison changes only the rank-transfer inference. Raw-margin versus
rank-margin results are:

| Policy | UMAP MCC | Scaffold MCC |
| --- | ---: | ---: |
| Original centered mean | .1005 / .0937 | .1765 / .1522 |
| Rank margin mean | .1162 / .0901 | .1470 / .1487 |

Mean rank-minus-original change is +.0060 UMAP, conditional interval
[-.0452,.0458], and -.0165 scaffold [-.0407,.0077]. The rank policy calls
28.3%/34.1% UMAP and 23.6%/36.6% scaffold. Aligning score scales neither establishes
a consistent MCC gain nor a stable call budget. Rank percentiles are not calibrated
class probabilities and are not assigned a Brier score.

## Connection to historical failures and the twelve-hour decision

The new evidence shows a concrete weakness of the new candidate: the chosen
decision boundary and its transfer vary substantially even when molecular scores
are strongly correlated. Extended fitting makes the score scale more extreme;
shorter fitting and rank normalization address that partially but do not establish
a robust decision improvement. This explains why counting 77 output positives or
quoting the network's earlier atomic MCC did not justify consuming an opportunity.

It does **not** identify the molecular cause of all historical Live failures.
The October 1 probe averaged ranks, so the new network's absolute probability-scale
problem cannot be asserted as that probe's direct cause. October 1 improved AP
without decision MCC improvement and exchanged identities at fixed 49 count.
Earlier AP-mass90 nevertheless failed after complete-pool and scaffold diagnostics
had supported it. Those are retained counterexamples to a universal AP, threshold-
scale, one-split or incomplete-pool explanation. Acquisition and scoring-mask
transfer remain unresolved, not solved by trying another architecture.

Actual-policy prior sensitivity at an illustrative .07 positive rate gives
original-policy UMAP MCC .0647/.0596 and scaffold .1149/.0960. The prior is unknown;
results at .03, .12 and the supplied training label fraction are also retained.
These are conditional label-shift scenarios, not Live predictions. Historical
whole-ensemble diagnostics retain the complete rate curves on exactly matching
IDs/folds; no convincing advantage over that whole ensemble is established.
Its evaluation analogy is itself imperfect and previously missed a Live reversal.

The submission would change the learner, decision policy, inference and call
count together, while the local evidence does not establish a repeatable gain.
It would provide a combined-pipeline observation, but little clean discrimination
among the failure hypotheses. Given the twelve-hour cost, retain the opportunity
instead of uploading this particular 77-call candidate. This neither adopts an
incumbent rollback nor fixes the next model to 49 calls. No replacement winner is
claimed from this controlled validation batch; the informative negatives remain.

## Reproducibility and retained records

The training/validation source is `scripts/validate_cyp_d6_bag_policy.py`, called
from the existing source-data evaluation through a recorded launcher. Input aliases
include canonical source tables, exact split assignments and historical complete-
pool comparisons. All substantive work used `pixi run cyp run`.

- `23751617`: original complete-policy validation and atomic controls.
- `74fc744f`: actual-policy prior, label-regime, threshold-scale and seed audit.
- `bdf341bd`: controlled 20-epoch validation.
- `0b8d2252`: 20-epoch audit, retained independently.
- `8e8ce151`: 20-versus-100 paired conditional intervals and audit.
- `17b6eddf`: controlled rank-transfer validation.
- `b0723c33`: rank-transfer prior/regime/seed audit.

All completed outputs are promoted to `artifacts/tracked/cyp-d6-zero-start-*`
with full run IDs, config and input checksums. Prepared payload
`20261002T023446563579Z-e264a4bf` remains immutable and unsent. Source checks pass;
infrastructure tests passed 73 with 8 skipped. Fit-boundary assertions and the
exact reproduction check verify the controls relevant to this experiment.
