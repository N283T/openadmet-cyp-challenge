# Comprehensive audit of weak CYP2D6 TDI predictions

The user's question was whether the earlier explanation was actually supported,
given previous implementation and interpretation mistakes. This report treats
that explanation as a hypothesis, not a premise. It covers classification TDI,
not CYP2D6 direct-inhibition regression. No production label, split, threshold
or submission was changed during this investigation.

## What the investigation changes

The statement that the saved D6 models have weak UMAP discrimination is supported.
The implied explanation that fewer labels account for the gap is insufficient.
The comparison with A4 conflated training size, different labelled chemical
populations and different observed-potency availability regimes. UMAP is also
substantially more structurally extrapolative than the test population.

No bulk label inversion, ID swap, endpoint mix-up or failed-to-fit model was found.
This does not establish a biological noise ceiling or absolve the modeling choices.
In particular, tree models fit their training rows almost perfectly and generalize
poorly on D6; the current weak result is not evidence that D6 cannot be modeled.

## 1. Verify the actual official comparison

Official leaderboard snapshots `738ff4c9`, `e4ea4a7a` and `64e5e66c`, included in
raw-provenance audit `e960da8e`, give:

| LB MCC | Current submission | Official TabICL | Leader on each board |
|---|---:|---:|---:|
| D6 | 0.176698 | 0.199337 | 0.430129 |
| A4 | 0.375036 | 0.450597 | 0.515446 |
| Macro | 0.275867 | 0.324967 | 0.450195 |

The official comparison is unfavorable, but D6's gap to TabICL is 0.022639, not
a massive discrepancy. The gap to the leading D6 submission is much larger.
The two endpoint leaders are different submissions; their scores are not a
real jointly submitted macro result. Published D6 bootstrap standard deviations
are 0.084810 (ours) and 0.090631 (TabICL). Separate standard deviations do not
provide a paired significance test. The Space reports bootstrap-mean metrics;
local pooled OOF MCC is a different estimator on different compounds.

Our old local baseline is TabPFN and is not the official TabICL baseline. Its
previous LB D6 MCC was 0.001644 and A4 0.430653. The new submission improves D6
and harms A4. It is wrong to describe that within-project improvement as having
beaten the official baseline.

## 2. Raw labels, identity and implementation checks

Audit `e960da8e` joins raw TRAIN_Emax and TRAIN_TDI records through canonical
aliases and checks the supplied booleans against the training labels. The 10,162
retained raw labelled occurrences resolve to 5,081 compound/endpoint labels:
1,497 D6 and 3,584 A4. There are zero conflicting labels, zero value mismatches,
zero unmapped retained source rows and no multi-name or multi-SMILES merges in
this classification-labelled cohort. Missing labels remain missing.

One historical exclusion matters for exact counts: raw Emax contains 1,498 D6
labels, but the platinum compound OCNT-2329004, `N.N.[Cl][Pt][Cl]`, was quarantined
by molecular standardization. Its D6 call is negative. The TDI table and canonical
classification training set contain 1,497. This is not evidence of a bulk label
problem, but saying every raw D6 label was retained would be inaccurate.

Follow-up `1edf91e1` distinguishes raw-label agreement from rule coverage. All
1,088 supplied positives reconstruct from the tested published potency rules,
but 22 D6 and 1,409 A4 supplied negatives are not reconstructable from that simple
rule and the available potency pairs. These include missing measurements and
18/159 complete-pair negatives outside the rule's defined branches. Zero
disagreements on covered rows must not be described as a complete independent
reproduction of every official call. The supplied calls are preserved; this
coverage limitation alone does not establish that they are wrong.

The previous independent audits checked saved LGBM booster predictions and inner
thresholds (`8aeb9a25`), TabPFN labels/fold masks and probabilities (`401ba682`),
joint fine-tuning masks/scalers/selected epochs (`7e38be86`), and all ensemble
probabilities, medians and binary outputs (`da4f1375`). The current diagnostic
`6559a18f` independently requires all 25 saved OOF methods to match canonical
compound/endpoint/fold/label keys. Probability columns take precedence over
binary columns; no orientation flip is applied. The submitted package `96483cf8`
also matched all 750 official names and original SMILES.

These checks cover identified failure modes; reproducing a prediction does not
prove that the chosen model architecture, capacity or auxiliary objective is good.

## 3. Weak OOF ordering is real, but not equivalent to no signal

Diagnostics `6559a18f` and `1eefaddd` evaluate 12 LGBM recipes, three TabPFN
recipes, two graph recipes, seven fixed means and the old baseline.

For the submitted D6 broad mean:

- Pooled AUROC is **0.566186**. Weighting only within-fold positive/negative
  comparisons gives **0.566514**: cross-fold score offsets do not explain it.
- Individual fold AUROCs are **0.6533, 0.5843, 0.5501, 0.4912, 0.5778**. Weakness
  spans multiple folds, while fold 3 is especially poor. The old model's fold 3
  AUROC is also only 0.4770.
- Conditional 95% bootstrap range for AUROC is **[0.5303, 0.5984]**. Within-fold
  label permutations preserve prevalence; observed ranking exceeds that null.
  Thus calling the model literally random would be too strong.
- Average precision is **0.2898**, against positive prevalence **0.2164**. The
  highest score decile has 34.9% positives. There is some enrichment, but limited.
- Brier score is **0.16951**, nearly the training-prevalence-only comparator
  **0.17007**. Log loss is **0.52626**, slightly worse than the comparator
  **0.52379**. The probabilities are not especially informative or well calibrated.

The D6 AUROC increase versus our old model is only **0.02333**, with paired
conditional bootstrap range **[-0.00658, 0.05305]**. The large LB MCC increase
therefore should not be interpreted as proof of a comparably large gain in
underlying ranking. For A4 the new-minus-old OOF AUROC change is **-0.01538**,
conditional range **[-0.02598, -0.00524]**. Old A4 should have been retained as a
serious alternative before submitting the replacement.

All uncertainty here is conditional on saved fits and already screened recipes.
Row resampling is not a refit-based or fully nested model-selection evaluation.

## 4. Thresholds matter, but cannot explain the whole OOF failure

The D6 selected inner cutoffs are 0.23065, 0.37324, 0.34458, 0.29465 and 0.26755.
They produce between one and 42 positives per outer fold. Such variability is
consistent with unstable threshold choice on small inner sets.

The observed inner-threshold OOF MCC is 0.09623. Even selecting the single best
threshold directly on these same outer OOF labels gives only **0.13219**. That
is an optimistic fitting diagnostic, not a deployable CV gain or a general
upper bound for other models. It shows that a better single monotone cutoff on
this particular score vector cannot create a strong OOF classifier.

The submission uses the median inner cutoff on fivefold mean probabilities.
This differs from one held-out model per OOF row. D6 mean probabilities are
0.1680 in OOF and 0.1729 on all test compounds; their standard deviations are
0.0602 and 0.0676. There is no gross average probability shift here. Alternative
majority voting would mark 46 rather than 43 test compounds positive and changes
19 individual calls. This is a diagnostic only; no replacement was deployed.

These are all-test statistics; the hidden confidently scored subset is unknown.
They cannot establish test calibration. Changing both the model and threshold
between submissions also prevents attributing the observed LB gain to either alone.

## 5. Controlled training-size experiments

Use balanced LGBM, Morgan + Mordred, identical 400-tree/15-leaf parameters and
the original UMAP assignments. The full controls reproduce stored predictions
to 1e-12. Subsampling changes only A4 outer-training rows, stratified by remaining
fold and label with three fixed seeds; the A4 validation set is unchanged.
This avoids confusing a smaller test sample with less training information.

| Training population | Pooled OOF AUROC |
|---|---:|
| D6, all labels | 0.533285 |
| A4, all labels | 0.764995 |
| A4, D6-matched training count, three seeds | 0.716742–0.742281 |
| A4, half the D6 training count, three seeds | 0.679709–0.711849 |

Reducing the count does hurt A4, but nowhere near the D6 result on that fixed
A4 validation population. UMAP-group bootstrap estimates the full-minus-matched
AUROC difference at 0.03518 [0.02349, 0.04577]. This is evidence against count
alone explaining the full gap; it does not isolate a causal endpoint effect.

Train/validation metrics also answer whether the trees learned anything. D6
mean training AUROC is **1.000**, training log loss **0.0634**; mean validation
AUROC is **0.5367**, log loss **0.5895**. A4's corresponding values are approximately
1.000/0.1363 and 0.7653/0.4659. D6's problem is generalization in this experiment,
not failure to optimize the training data.

A deliberately stronger-regularization control (`465afc7e`: 100 trees, seven
leaves, minimum child 80, L2 10) gives D6 OOF AUROC **0.530677**, compared with
the original **0.533285**. A shuffled-training-label control with the original
capacity still reaches training AUROC approximately 1.0, but true OOF AUROC
**0.509864**. Therefore the near-perfect training score demonstrates memorization
capacity, not learning a transferable TDI mechanism. The tested regularization
change does not recover generalization. One label permutation per fold is a
sanity control, not a statistical significance distribution or exhaustive tuning.

## 6. The endpoint label populations are materially different

The overall positive rates look similar, but mask a major composition difference:

| Supplied labels | D6 | A4 |
|---|---:|---:|
| Both observed direct and TDI pIC50 available | 1,493 | 2,334 |
| Incomplete potency pair | 4 | 1,250 |
| Positives in incomplete-pair rows | 0 | 0 |

Absent canonical potency means no usable value in this table; it does not mean
the assay was never performed. The extra A4 negatives constitute an easier
regime for the present classifier. Restricting evaluation of the saved selected
A4 mean to both-measured rows lowers AUROC from **0.7779 to 0.6950**. D6's
both-measured AUROC remains **0.5665**.

The additional training control restricts A4 training AND evaluation to its
2,334 both-measured compounds. Balanced FP+Mordred AUROC is **0.708457** with
all such training labels, and **0.683303–0.688511** after matching D6's training
size. A gap remains even after these two controls, but these are still different
chemical populations.

Only **259 compounds** have labels for both endpoints. Applying the saved
full-population ensembles to this same subset gives D6/A4 AUROC **0.6131/0.6430**.
Training both endpoints using only this shared population gives **0.5634/0.4769**;
the apparent reversal has a broad paired group-bootstrap interval crossing zero.
The shared subset also has 19.7% D6 versus 40.9% A4 positives. Therefore the data
do not support an invariant claim that D6 is intrinsically harder in every cohort.

## 7. Structural coverage and label disagreement

Using exact Morgan radius-2/2048-bit Tanimoto, D6's nearest labelled training
neighbor has median similarity **0.2985** in UMAP OOF, versus **0.4711** for test
against all training labels. This is not just caused by having more candidate
neighbors at test inference: with identical per-fold candidate pools, test
medians remain **0.380–0.409**, versus OOF **0.292–0.317**. Test fractions with
similarity at least 0.5 are **34–42%**, versus OOF **0.35–1.83%**. A4 shows the
same direction of difference.

Thus this UMAP validation setup evaluates much stronger structural extrapolation
than is typical for test. It is plausible that this contributes to OOF/LB
differences; similarity does not reveal the hidden labels or quantify the score
effect. The split has not been silently changed.

Nearest-neighbor label disagreement is **31.26% for D6 and 31.42% for A4**, versus
within-fold prevalence-preserving null means 33.90% and 33.59%. Most neighbors are
not highly similar, and very close OOF pairs are too sparse for strong claims.
This analysis does **not** support calling D6 labels uniquely noisy. Measurement
noise, activity cliffs and inadequate fingerprint relevance cannot be separated
with these disagreement counts alone.

Errors are not confined to weak or ambiguous observed shifts: 244 of 276 positive
D6 examples with strictly positive conservative shift bounds are missed at the
saved cutoff. Observed measurements are used only for retrospective strata, not
as inference features or replacement labels.

## What is supported, and what remains open

Supported: weak structural generalization of these D6 models under original UMAP;
limited but nonzero ranking information; unstable cutoff selection; different
labelled cohorts; an A4 negative regime largely absent from D6; a structural
OOF/test mismatch; and a training-size effect insufficient to explain all of the gap.

Not established: a corrupted dataset, uniquely noisy D6 biology, an unavoidable
D6 performance ceiling, a threshold-only cause, or failure of external data in
general. The tested neural objectives and fixed hyperparameters are a limited
set of methods, not an exhaustive search of useful TDI representations.

The external encoders `a3ea2ec3` and `2c0a0953` were supervised with ChEMBL
potency/activity, PubChem activity and Tox21 AC50/activity endpoints. They were
not trained on matched direct-versus-preincubated TDI shifts. Their failure here
does not test every form of external TDI supervision, nor establish that such
supervision would necessarily help.

Next useful comparisons should explicitly control the labelled chemical/assay
populations, retain the original UMAP as a reference, and evaluate whether a
test-coverage-matched diagnostic split changes model rankings. Representation or
auxiliary improvements should target TDI-specific distinctions, with valid
cross-fitting; simply adding more direct-inhibition encoders is not a demonstrated
solution. These are proposed next experiments, not completed improvements.

## Recorded evidence

| Run suffix | Work |
|---|---|
| `e960da8e` | Raw-label reconciliation, rule checks and official scores |
| `1edf91e1` | Rule coverage versus disagreement and external-objective audit |
| `6559a18f` | 25-model ranking/calibration/fold/conditional uncertainty audit |
| `1eefaddd` | Paired AUC changes, cutoff sensitivity and deployment reconstruction |
| `5bf02051` | 40 new size/shared-cohort fits and 10 reproduced full controls |
| `f551edca` | UMAP-group uncertainty for size/shared comparisons |
| `d14ded5f` | 20 both-assay A4 size-control fits |
| `2a487983` | Final restricted-size and train/validation summary |
| `4d006374` | Structural neighbors, label routes, errors and permutation controls |
| `7f106e62` | Identical-pool OOF/test structural comparison |
| `465afc7e` | Ten regularization and shuffled-training-label sanity fits |

See [neighborhood and label routes](2026-09-12-classification-d6-neighborhood-routes.md),
[controlled learning comparisons](2026-09-12-tdi-learning-size.md) and
[training sanity controls](2026-09-12-classification-d6-training-sanity.md), plus
[raw provenance and official evaluation](2026-09-12-classification-label-provenance.md).
Run folders retain complete configurations, declared input snapshots, selected
training compound IDs, trained models and keyed predictions. All analysis and
training ran through `pixi run cyp run`; raw snapshots and prior artifacts remain intact.
Raw-provenance attempt `dc511012` failed on the initial overly strict assumption
that every raw row has a canonical ID. It remains recorded; `e960da8e` explicitly
accounts for the pre-existing quarantined compound. No reported result relies
on the failed attempt. Final Ruff checks passed; infrastructure tests completed
with 37 passed and seven skipped. An independent review of the combined report
checked the learning-study values and statistical qualifications against the
underlying artifacts.
