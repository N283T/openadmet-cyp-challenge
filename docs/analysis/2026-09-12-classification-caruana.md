# Classification ensemble selection and compact alternatives

The user asked us to carry out the previously authorized Caruana comparison,
including reducing an unnecessarily large classification ensemble. The aim is
to test selection and weighting on held-out predictions, not assume they help.

## Selection protocol

Every outer UMAP fold uses the designated next fold as inner validation. Its
candidate predictions come from fits excluding both the outer and inner folds.
Weights, model identities, ensemble length and thresholds are selected using
only these inner predictions and labels. The untouched outer fold evaluates
that entire selection procedure. The same original label/split keys are checked
across all candidates. This is not selection on other folds' ordinary OOF
predictions from models that already trained on the outer labels.

The compared procedures are:

- Best single model, exhaustive equal-weight pair, exhaustive equal-weight
  triple, each selected by inner MCC with an inner-selected threshold.
- Caruana with replacement, optimizing MCC, AUROC or negative Brier score,
  up to 50 additions. Retain the earliest strictly best prefix. MCC scoring
  searches unique inner scores plus 0.5, using `>=`, with ties closest to 0.5
  and then the larger cutoff. Stable alphabetical model order resolves equal
  candidate scores.
- MCC Caruana restricted to at most three distinct model types.
- In the expanded pool: fixed historical TabPFN control and bagged sparse
  Caruana. The latter runs sparse Caruana on 20 class-stratified bootstrap
  samples, averages the weights, keeps the largest three and normalizes them.
  Its final cutoff is fitted on the original inner set. This is an attempted
  stabilization, not a claim that bagging must improve performance.

For final test artifacts, each recipe fits weights and a cutoff to all OOF
labels, then applies them to the average of each candidate's five test-fold
probabilities. These fitting scores are stored separately and never called CV.
The outer-evaluation weights and fold-test probabilities are also preserved.
No all-label base-model refit or new remote submission is implicit in this work.

## Initial 17-candidate pool

Run `d27648b2` uses twelve LGBM, three frozen TabPFN and two fine-tuned graph
classifiers with previously audited inner/outer predictions. Existing fixed
references are the submitted D6 eleven-type mean and A4 balanced six-LGBM mean.

| Selection procedure | D6 outer MCC | A4 outer MCC |
|---|---:|---:|
| Existing fixed reference | 0.096228 | 0.338112 |
| Best single | 0.042403 | 0.293298 |
| Best equal-weight pair | 0.051188 | 0.318775 |
| Best equal-weight triple | 0.050261 | 0.293480 |
| MCC Caruana | 0.047697 | 0.317735 |
| MCC Caruana, at most three types | 0.061665 | 0.331549 |
| AUROC Caruana | 0.070353 | 0.310563 |
| Brier Caruana | 0.038395 | 0.307723 |

The initial pool does not justify replacing the fixed references. In contrast,
full-OOF-fitted MCC Caruana scores are 0.187301/0.391490: reporting those as
validation improvement would be misleading. Small inner labelled sets allow
selection and threshold choices to fit noise; the outer results test the actual
consequence. Summary run `30f9ee8f` preserves paired UMAP-group intervals and
weight complexity. Its intervals are conditional on fitted and screened models.

Independent audit `f04fde09` replays every greedy step, exhaustive subset,
best-prefix choice and cutoff, checks the actual canonical inner fold, and
reproduces 35,567 OOF, 52,500 fold-test and 10,500 final test rows. Initial audit
attempt `e00ecc2e` failed a nested filename lookup ambiguity, fixed by exact
declared-name priority; it supplied no reported result and remains recorded.

## Historical component inclusion

The historical stronger A4 model is PCA256 of frozen CheMeleon plus four nested
predicted log2FC values plus Mordred, followed by TabPFN 2.6 with eight estimators.
Its old outer predictions and test-fold probabilities are retained unchanged.

Existing pair-excluded log2FC models are insufficient for its new inner fits:
readout-training fold `r` needs an auxiliary encoder excluding `{outer, inner, r}`.
Reusing the old stacked matrix would expose held-out inner auxiliary labels.
The exact reconstruction therefore trains ten triple-excluded auxiliary models
and ten inner classifiers, with train-only PCA and Mordred preprocessing.
These auxiliary fits are evaluation dependencies, not extra models silently
added to the final deployed candidate. The historical deployment retains its
own saved pair-excluded auxiliary models.

Reconstruction `130ac6af` completed all ten auxiliary and ten readout fits;
independent provenance audit `68fe5a1c` passed. See the
[historical reconstruction report](2026-09-12-historical-stack-inner.md).

## Completed 18-candidate comparison

Run `8e3be966` adds this historical candidate, with summary `d19c78b7` and
independent selection replay `22321032`. All nine procedures completed.

| Procedure | D6 outer MCC | A4 outer MCC |
|---|---:|---:|
| Existing fixed submitted reference | 0.096228 | 0.338112 |
| Fixed historical TabPFN, honest inner cutoff | 0.026018 | **0.365149** |
| Inner-selected single | 0.042403 | 0.335881 |
| Inner-selected equal pair | 0.051188 | 0.316770 |
| Inner-selected equal triple | 0.050261 | 0.329967 |
| MCC Caruana | 0.071475 | 0.337668 |
| MCC Caruana, maximum three types | 0.061777 | 0.356736 |
| AUROC Caruana | 0.074920 | 0.342891 |
| Brier Caruana | 0.038395 | 0.344491 |
| Bagged Caruana, maximum three types | 0.084431 | 0.332283 |

The strongest tested A4 Caruana procedure is the three-type cap, improving the
submitted fixed reference by 0.018625; its conditional paired group-bootstrap
interval is [-0.005316, 0.045415]. It still does not exceed fixed historical
TabPFN. Its final full-OOF fit happens to retain only two types, equally weighted:
historical stack and balanced external-classification-embedding + Mordred LGBM.
The displayed outer score evaluates the selection procedure, not this exact
fixed pair chosen on all OOF labels.

D6's best tested Caruana procedure is bagged sparse selection, still below the
existing fixed reference. Its final weights are external-mixed LGBM 0.390515,
Morgan balanced LGBM 0.322752 and external-mixed PCA256 TabPFN 0.286733. It is
not adopted merely because it has fewer models.

Unrestricted MCC Caruana's full-OOF D6 fitting MCC is **0.205508**, versus honest
outer MCC **0.071475**. This sharp distinction is why full-OOF fit numbers must
not be reported as an improvement. Model selection from only one designated
inner fold is valid with respect to the outer labels, but its small sample is
a limitation. This work does not claim that all forms of ensemble selection
or richer inner cross-validation would fail.

The independent expanded audit reconstructs all 108 selections, including all
20 resampled selections in each bagged fit, and verifies 45,729 OOF, 67,500
fold-test and 13,500 deployment rows plus both outer and fitting-only metrics.
See [selection audit](2026-09-12-classification-caruana-audit.md).

## Compact fixed-model alternative and prepared files

The earlier fixed Morgan balanced LGBM with cutoff 0.5 has D6 OOF MCC
**0.093546**, compared with **0.096228** for the eleven-type reference. Packaging
run `effedde5` verifies that exact saved prediction comparison and a conditional
group-bootstrap difference interval [-0.07561, 0.08850]. This does not prove
equivalence; the fixed recipe was identified in the prior screen, unlike the
inner-selected-single procedure above. It supplies a practical one-type
comparison, with five readout models instead of 55.

Two complete 750-row files were prepared and dry-run validated:

| File | D6 | A4 | Positive calls D6/A4 |
|---|---|---|---:|
| `retain_d6_restore_a4` | Exact current submitted vector | Exact previous submitted vector | 43 / 249 |
| `compact_morgan_restore_a4` | Morgan balanced fivefold mean, cutoff 0.5 | Exact previous submitted vector | 117 / 249 |

Both preserve official IDs/SMILES and the previous A4 boolean column exactly.
The compact version changes many test calls; similar OOF MCC does not imply
identical LB behavior. The A4 historical stack also has shared auxiliary model
dependencies, so five readouts does not mean five total neural assets.

Candidate `e945c8b7` / prepared `0281edac` is the exact restoration option.
Candidate `c9827313` / prepared `8955b0f0` is the compact comparison. Their immutable
prediction sources and CSVs are under run `20260912T133350390925Z-effedde5`.
Neither was uploaded. All Caruana recipes also retain their own full-OOF-fitted
test probabilities, weights and cuts, even when they were not adopted.

Decision before the TabICL expansion below: the comparison includes the strongest historical
candidate. Caruana does not currently justify replacing D6's fixed reference or
A4's historical single stack. Keep the compact D6 option as a concrete next
comparison, rather than treating an eleven-type ensemble as necessary.

## TabICL-expanded 24-candidate comparison

Run `97df4c71` adds three frozen-view TabICL candidates, matched Jeremy600
TabICL8 and TabPFN8, and the historical-stack TabICL replacement. The original
18 candidates, nine procedures and exact
`history::historical_stack_mordred` control remain unchanged. Test probabilities,
outer-weighted fold predictions, deployment weights/cutoffs and all separate
full-OOF fitting values are saved.

| Procedure | D6 outer MCC | A4 outer MCC |
|---|---:|---:|
| Existing fixed submitted reference | **0.096228** | 0.338112 |
| Fixed historical TabPFN | 0.026018 | 0.365149 |
| Inner-selected single | 0.042403 | 0.305745 |
| Inner-selected equal pair | 0.044502 | 0.348422 |
| Inner-selected equal triple | 0.078507 | 0.330445 |
| MCC Caruana | 0.088698 | 0.346761 |
| MCC Caruana, maximum three types | 0.053293 | 0.349099 |
| AUROC Caruana | 0.092594 | **0.387067** |
| Brier Caruana | 0.030865 | 0.335858 |
| Bagged Caruana, maximum three types | 0.050436 | 0.331489 |

The new candidates do not improve D6 beyond its fixed submitted reference.
A4 AUROC-Caruana improves the submitted fixed reference by 0.048956 and the old
historical stack by 0.021918. Summary run `8889bb1f` gives a conditional paired
UMAP-group interval of [0.019969, 0.077007] against the submitted A4 reference.
The all-OOF deployment fit uses six base types (30 fold readouts), cutoff
0.233569 and 361 blind-test positives. Its full-OOF MCC 0.417948 is a fitting
value, not the honest outer score.

Independent audit `d8ca5f7b` reconstructs all 108 selections and metrics for
24 candidates, checks canonical inner/outer folds and labels, and reproduces
all saved OOF, fold-test and deployment predictions. The practical target-wise
candidate retains the current D6 vector and replaces only A4 with this ensemble;
packaging run `d30a0365`, candidate `f0deac96` and prepared submission `82d53182`
passed validation and dry run. No upload was made.
