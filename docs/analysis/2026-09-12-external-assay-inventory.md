# Assay-aware external CYP inventory and overlap diagnostics

The combined external-only dataset contains **35,320 standardized compounds,
156,372 source observations, 67,824 exact continuous labels and 125,063 binary
labels**. It combines ChEMBL 36, PubChem AID1851 and five Tox21 P450-Glo assays.
All challenge connectivity blocks are excluded from fitting data. Source-specific
assay identities and label semantics remain separate.

## Immutable runs and usable assets

- Source schema inspection: `20260912T050756113457Z-89b0a428`.
- Final ChEMBL/PubChem curation: `20260912T051342149959Z-2b7ef7dd`.
- Tox21 source curation: `20260912T051221903183Z-6113ed82`.
- **Combined dataset: `20260912T051513413487Z-ac9df6b5`**.
  Use `external-compounds.parquet` and `observations.parquet` in its outputs.
  `inventory.csv`, `binary-conflicts.parquet` and `summary.json` document coverage.
- Curation label/identity audit: `20260912T051150285630Z-5bc5535c`.
- Independent loss audit: `20260912T051644620106Z-fca10454`.
- **External/challenge overlap diagnostic: `20260912T051816691443Z-ba507a53`**.

Configuration: `configs/external-assay/curation-v3.json` and
`configs/external-assay/merged-tox21-v1.json`. Implementations:
`scripts/external_assay_prepare.py`, `scripts/external_assay_merge_tox21.py`,
`scripts/external_assay_overlap.py`.

## Label and identity policy

ChEMBL selection requires a human single-protein target, human or unspecified
assay organism, confidence >=8, A/B/F assay type, no validity/duplicate flag and
no explicit TDI/time-dependent/mechanism-based/irreversible description. Explicit
preincubation is flagged and retained; it is not automatically called TDI. Unlike
the old target-organism-only filter, explicitly nonhuman assay organisms are
excluded. All ChEMBL PubChem imports are excluded to avoid repeating AID1851.

IC50 and Ki remain different measurement types. All accepted nM values convert
as `9 - log10(nM)`. The inequality reverses under this transform. Censored values
have null `value` and preserve `lower_bound`/`upper_bound`; they are not exact
potencies. A ChEMBL binary label uses p-potency >=5 only when exact values or
bounds establish the class unambiguously. Original inactive comments can supply
binary-only labels when no numeric potency exists.

PubChem AID1851 continuous labels use `-Fit_LogAC50`, whose source definition is
log10 molar AC50. Inactive observations never receive potency values. Original
Active/Inactive calls supply binary labels; Inconclusive is missing binary, and
its finite fitted potency remains identifiable by hitcall/curve metadata.
Unclipped Max_Response, fit R2, curve class and source QC remain separate columns.

Tox21 binary labels follow original overall Active/Inactive calls. Potency uses
median replicate pAC50 only for overall Active, inhibitor phenotype and curve
classes -1.1, -1.2, -2.1 or -2.2. Source run retains original SID rows and replicate
metadata. The merge copies the supplied aggregate target table for provenance.

The merge masks **1,047 contradictory compound/assay binary keys**, including
fractional means of disagreeing Tox21 SIDs. `original_binary_label` and
`binary_conflict` preserve the decision. No contradictory class is silently
converted to inactive. Exact compound identities are deduplicated across sources,
while observations and assays are retained separately.

## Coverage before Tox21 addition

The ChEMBL/PubChem base contains 28,993 compounds, 116,360 observations,
57,735 exact continuous labels and 17,014 censored measurements. It spans 8,745
source/assay/CYP/type keys. There are 1,963 assays with at least two distinct
compounds and continuous values; 769 contain at least five compounds.
These are inventory eligibility counts, before the trainer's minimum-size and
training/validation criteria. Tox21 adds 6,327 unique compounds to this base.

## Agreement with challenge measurements

The diagnostic uses exact standardized identities as the primary comparison.
It finds **195 compounds and 458 compound/assay/CYP pairs**, across 87 assays.
Only 15 assays have >=3 shared compounds and 9 have >=10. Connectivity matching
adds 175 nonidentical pairs, which are saved separately rather than treating
stereoisomers as identical calibration observations.

`offset` below means **challenge pIC50 minus external p-potency**. These statistics
are descriptive and have not been applied to any OOF prediction.

| External assay | CYP | Shared compounds | Spearman | Median offset |
|---|---|---:|---:|---:|
| PubChem1851 | 1A2 | 32 | 0.7715 | -0.7041 |
| PubChem1851 | 2C9 | 8 | 0.1687 | -0.9654 |
| PubChem1851 | 2D6 | 46 | 0.7634 | -0.3562 |
| PubChem1851 | 3A4 | 10 | 0.5461 | -0.8770 |
| Tox21 | 1A2 | 57 | 0.6243 | +0.1624 |
| Tox21 | 2C9 | 22 | 0.4589 | -0.0865 |
| Tox21 | 2D6 | 114 | 0.8275 | -0.5200 |
| Tox21 | 3A4 | 24 | 0.5842 | -0.8336 |

For 2D6, both panels preserve substantial ordering despite offsets. For 1A2,
Tox21 and AID1851 offsets differ markedly; one global external correction would
not match both sources. Assay-specific or source-specific training therefore has
empirical support. ChEMBL source-level statistics repeat compounds measured in
multiple assays; they are not independent-compound correlations. Most individual
ChEMBL assay overlaps are too small to estimate an unrestricted slope reliably.

A future calibration experiment can estimate source/assay offsets or shrinkage
within each challenge training fold. Current diagnostics use available challenge
labels for description only and do not establish calibrated OOF performance.

## Trainer review and independent mathematical checks

The reviewed trainer separates continuous/ranking potency outputs from binary
logits, uses connectivity-group validation shared across sources, and fits scales
on external training values only. Ranking comparisons remain inside an assay and
require nonoverlapping potency intervals. Source-balanced training sampling is
implemented. Replicate intervals intersect (null bounds represent an unbounded
side), inconsistent intersections are masked, and exact replicates use medians.

Independent CPU checks validate 128 generated ordered pairs, ranking gradient
direction, and one-sided interval gradients including unbounded sides. All pass.
Review identified checkpoint averaging across assays as inconsistent with
source-balanced training. The trainer owner corrected validation to balance
sources before the full program launched; its `validation_by_source` logs expose
the components. These checks concern implementation, not predictive quality.

## Preserved revisions and failures

- Curation v1 `20260912T051008034379Z-e692eb5b` was superseded for strict
  threshold-bound binary handling and regex review.
- V2 `20260912T051120437541Z-e590c789` still contained a double-escaped TDI word
  boundary. Recorded regex audit `20260912T051309231712Z-8591d209` failed on
  literal `TDI`, exposing that the initial edit had not corrected serialization.
- Attempt `20260912T051322189350Z-7ea5ef84` failed because the v3 configuration
  file had not been created; the shell's `python` was unavailable.
- V3 regex audit `20260912T051341612033Z-4732d338` passes literal TDI tests. It
  identifies 31 explicit-TDI source rows. The corrected exclusion does not change
  final accepted counts because affected rows also fail other selection criteria.
- All earlier outputs remain intact; final merged inputs use corrected v3.

## Subsequent NVS addition (separate frozen corpus)

Run **`20260912T052252244003Z-ab6072d6`** adds NVS_ADME assays to the
three-source corpus without changing the inputs of the initial four-objective
comparison. Its output schema is the same `external-compounds.parquet` /
`observations.parquet` contract.

The resulting corpus has **35,356 compounds, 158,027 observations, 68,332
continuous labels and 126,208 binary labels**. NVS contributes 590 compounds,
of which 554 already occur in the prior corpus: only **36 are new identities**,
while 1,655 new source measurements contribute 508 continuous and 1,145 usable
binary labels. Two additional contradictory binary keys are masked, bringing
the total to 1,049. This is primarily assay enrichment, not a large chemical
coverage expansion.

The source is EPA invitrodb v3.3 (2020), with CAS identities mapped to the official
2018-10-17 DSSTox structure snapshot; source curation run
`20260912T052128343558Z-128af2eb` records its provenance. Binary labels require
representative `gsid_rep=1`, `hitc` in {0,1} and micromolar concentration units.
Potency additionally requires active, finite, unflagged observations and converts
`modl_ga = log10(uM)` as **`pAC50 = 6 - modl_ga`**. NVS assay IDs remain distinct
from PubChem, Tox21 and ChEMBL. Original assay names, tested concentration units,
response units and curated source records are included in merged run outputs.
All challenge connectivity exclusions and prior binary masks are verified again.

## Small-assay normalization audit

Recorded CPU audit **`20260912T052617517710Z-c45b7793`** inspects task metadata
and training membership of the continuous run
`20260912T052347954512Z-91ba5769`. It reconstructs the trainer's replicate policy
and asserts that every training mean and standard deviation matches saved task
metadata. Outputs include `assay-scale-audit.parquet`, `source-scale-summary.csv`,
`train-only-pooled-scales.csv` and `high-leverage-assays.csv`.

| Source | Continuous assays | Median exact training n | Median sigma | Sigma at 0.1 floor | Mean-predictor raw validation interval MSE | Normalized MSE |
|---|---:|---:|---:|---:|---:|---:|
| ChEMBL36 | 996 | 7 | 0.397 | 97 | 0.446 | 3.419 |
| PubChem1851 | 5 | 8,172 | 0.700 | 0 | 0.423 | 0.938 |
| Tox21 | 5 | 1,916 | 0.577 | 0 | 0.355 | 1.021 |

The predictor in this table is each assay's **training mean**, not a neural model.
The comparison isolates the scale convention: ChEMBL physical-unit errors are
similar to the larger panels, while small fitted sigmas inflate normalized errors.
Among ChEMBL tasks, 405 have <=5 exact training values and 263 have sigma <=0.25.
The 90th percentile inverse variance is 96.4. One task has two exact training
values spanning only 0.173 p-units; its sigma reaches the 0.1 floor and a raw
validation squared error of 1.41 becomes a normalized error of 141.1.

A prospective **optional** scaling variant should retain assay-specific training
means but regularize variance using training-only source × CYP × measurement-type
pools. A within-assay residual pool avoids interpreting assay mean shifts as
within-assay variability. One declared candidate is

`sigma² = (n_exact * assay_variance + 20 * pooled_within_assay_variance) / (n_exact + 20)`,

followed by a sigma floor of 0.3. The audit computes this candidate using only
training exact values; validation values never define scales. For the same mean
predictor, ChEMBL normalized MSE becomes 1.615 under that convention. **This is a
change in loss units, not an accuracy improvement.** The appropriate comparison
is downstream OOF quality after an independently recorded training run. Existing
four-objective comparison configurations and running models are unchanged.

The optional `scaling_policy="pooled_shrinkage"` is now implemented in
`scripts/train_external_assay.py`, with helper `external_assay_scaling.py`.
The default remains `"assay"` and preserves its numeric scaling behavior.
`scaling_pseudocount=20` and `scaling_minimum_std=0.3` are explicit effective
configuration fields. Pools are calculated once, before task construction, from
training exact deduplicated labels. Per-assay means remain unchanged. Group,
prior variance, training count, fallback status and shrinkage parameters are
saved in task/checkpoint metadata; the complete pool table is a run output.

Recorded verification **`20260912T052936486765Z-7e14a6cc`** perturbs every
held-out lower/upper label as `-1000*x + 100000` and confirms the 20 variance
pools are exactly unchanged. Perturbing a training assay changes its pool,
checking the test is not vacuous. The helper's pools also match the prior
independent audit across all 1,006 eligible continuous assays, and direct
formula/floor cases pass. No current experiment configurations were modified.
