# Controlled TDI learning-size and population comparisons

Reducing3A4 training to the2D6 sample count lowers its AUROC, but does not reproduce
the weak2D6 result. Restricting3A4 to compounds with finite measurements in both
assay arms reduces the advantage further, yet a gap remains. On the very small
population labelled for both targets, the numerical ordering reverses. These
results support checking the label/chemical population and generalization gap,
rather than attributing the whole difference to sample count or an intrinsically
harder endpoint.

## Controlled design and provenance

The fixed learner is the balanced Morgan2048(radius2)+Mordred LGBM recipe from
`76facc71`:400 trees,15 leaves, minimum child20, learning rate0.03, L2=2, seed42+fold.
The original five UMAP folds remain unchanged. Ten saved full-data fold boosters
serve as exact controls; reconstructed validation probabilities match their saved
OOF values within1e-12. Their train probabilities are also evaluated.

Run `20260912T125841153686Z-5bf02051` fits40 new models:

- For each3A4 outer fold, its full validation population remains fixed. Training
  is reduced to that fold's2D6 training count, separately at seeds17/42/73.
  Sampling is proportional within training-fold×class strata, with exact total
  allocation using largest remainders. A nested half-size sample is drawn from
  each matched subset. No held-out labels select models or sample seeds.
- Both endpoints are separately trained and evaluated on the exact same259
  compounds with both target labels, preserving the original folds.

Run `20260912T130141371276Z-d14ded5f` adds20 models on3A4's2,334-compound
both-assay-measured cohort from route audit `4d006374`:full training and three
D6-sized subsets, with the same restricted validation population held fixed.
“Both measured” means finite direct and TDI potency in the available records;
missing potency is not interpreted as an assay never having been performed.

Every selected training ID, label, source fold and outer fold is saved. Audits
`f551edca` and `2a487983` verify fold exclusion, exact sample sizes, nested half
membership, fixed validation identities and identical shared-population IDs.

## Fixed-population validation results

Values below are pooled OOF metrics. Seed ranges describe three separate fitted
models, not an ensemble of their predictions.

| Target / validation population | Mean training rows per fold | AUROC | AP |
|---|---:|---:|---:|
| D6 full1,497 | 1,197.6 | 0.5333 | 0.2612 |
| 3A4 full3,584, full training | 2,867.2 | 0.7650 | 0.4464 |
| Same3A4 validation, D6-sized training | 1,197.6 | 0.7167–0.7423 | 0.3826–0.4219 |
| Same3A4 validation, half-D6-sized training | 598.4 | 0.6797–0.7118 | 0.3500–0.3910 |
| 3A4 both-assay2,334, full training | 1,867.2 | 0.7085 | 0.5257 |
| Same both-assay validation, D6-sized training | 1,197.6 | 0.6833–0.6885 | 0.4848–0.4941 |

The paired500-replicate bootstrap resamples UMAP chemical groups separately
within outer folds, using identical draws on both sides. Seed summaries average
metric differences across seeds, not probabilities. This quantifies uncertainty
conditional on these fitted models and folds; it does not refit each bootstrap.

| Fixed-validation contrast | AUROC loss from reducing training | 95% group-bootstrap interval |
|---|---:|---:|
| Full3A4 to D6 size | 0.03518 | 0.02349–0.04577 |
| Full3A4 to half D6 size | 0.06937 | 0.05458–0.08492 |
| D6-sized3A4 to half size | 0.03419 | 0.02340–0.04496 |
| Both-assay3A4 full to D6 size | 0.02221 | 0.01383–0.03029 |

The fixed learner benefits from additional data. However, even D6-sized3A4
training on the both-assay population reaches approximately0.686 AUROC,
versus0.533 for fullD6. These comparisons do not isolate all cross-target
differences: labels, compounds, chemistry and assay populations still differ.

## Training fit versus generalization

These are arithmetic means of fold-specific metrics, so they differ slightly
from the pooled OOF table above. Log loss uses the actual class labels without
reweighting the evaluation population.

| Population | Train AUROC | Validation AUROC | Train log loss | Validation log loss |
|---|---:|---:|---:|---:|
| D6 full | 1.0000 | 0.5367 | 0.0634 | 0.5895 |
| 3A4 full | 1.0000 | 0.7653 | 0.1363 | 0.4659 |
| 3A4 both-assay | 1.0000 | 0.7138 | 0.1319 | 0.5976 |

The model fits the D6 training labels extremely well. The weak validation result
therefore reflects poor generalization under this representation/model/split,
not an inability to optimize the training objective. This pattern alone does not
identify whether weak transferable signal, label construction, noise or chemistry
coverage is responsible.

## Shared259-compound comparison

Both endpoints use the same approximately207 training and52 validation compounds
per fold. D6 has51 positives(19.69%);3A4 has106(40.93%). Pooled AUROC is0.5634
forD6 and0.4769 for3A4. The paired3A4-minus-D6 AUROC difference is−0.0865,
with95% group-bootstrap interval−0.1875 to+0.0067: the observed reversal is
uncertain. AP values0.2368 and0.4022 have substantially different prevalence
baselines and should not be read as a direct superiority comparison.

This tiny selected cohort does not support a universal claim that D6 is harder
than3A4. Both shared-cohort training AUROCs are1.0, while validation remains weak.

## Artifacts

- `5bf02051`: original40 fits, ten exact controls, `selected-training-ids.parquet`,
  train/validation `predictions.parquet`, fold and pooled metrics.
- `f551edca`: original membership audit, paired group-bootstrap intervals and
  all500 replicate differences, mean fold metrics.
- `d14ded5f`:20 both-assay restricted fits and exact selected IDs/predictions.
- `2a487983`: all fold train/validation AUROC, AP and log loss; pooled prevalence;
  both-assay membership audit and paired bootstrap.

Scripts: `tdi_learning_size.py`, `tdi_learning_size_launch.py`,
`tdi_learning_size_summary.py`, `tdi_learning_final.py`. All experiments are
recorded classification runs with declared inputs. No submissions were made.
