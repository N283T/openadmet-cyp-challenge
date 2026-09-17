# Intrinsic CYP task design for the project-owned Chemprop encoder

## Question

The previous encoder trained 67 assay heads through 131 sampled objectives:
continuous interval regression, within-assay ranking and binary classification.
This screen asks whether that broad mixed objective improves Chemprop's own
held-out external-task accuracy. Challenge-label readouts are deliberately out
of scope.

All comparisons use connectivity-disjoint 10% validation splits, standardized
interval MSE per assay and the macro mean over the same 14 evaluable CYP
continuous assays. Lower is better. The external corpus excludes challenge
connectivity before model training.

## Objective and sampling screen

| Training design | Run | Heads | CYP task-macro MSE | Observation-weighted MSE |
|---|---|---:|---:|---:|
| Broad continuous | `17117688` | 54 | .57607 | .49485 |
| CYP continuous, source-uniform | `311ac3d5` | 15 | .52909 | .40896 |
| CYP continuous + weak ranking | `9605a976` | 15 | .53206 | .42312 |
| CYP continuous, task-uniform | `3a64d8b3` | 14 | .52883 | .41987 |
| CYP continuous, count^0.25 | `a583fc38` | 14 | .51443 | .40914 |
| CYP continuous, count^0.50 | `be0633d5` | 14 | .51157 | .39809 |
| CYP continuous, count^0.75 | `5a86f515` | 14 | .51732 | .39411 |
| CYP continuous, count^1.00 | `12f70712` | 14 | **.51065** | **.38565** |

The 15th CYP assay has only eight continuous observations and no useful held-out
support. Raising the eligibility threshold from five to 50 removes it. Sampling
assays proportional to their training observations performs best on both macro
and observation-weighted metrics. Adding ranking does not improve continuous
task accuracy.

## Head sharing and architecture

Run `74106180` retains 14 assay-specific normalized losses but shares one head
for matching CYP isoform and measurement type across sources. This reduces 14
heads to nine but worsens macro MSE to `.51485`; cross-assay head sharing is not
retained.

Architecture run `6543c9c0` fixes the selected 14-task continuous objective and
observation-proportional sampling. Its results are:

| Architecture | Task-macro MSE |
|---|---:|
| width600 / depth3 / FFN3 | .50669 |
| width600 / depth2 / FFN3 | .49983 |
| width600 / depth4 / FFN3 | .50576 |
| width448 / depth3 / FFN3 | .50401 |
| width768 / depth3 / FFN3 | .51230 |
| width600 / depth3 / FFN2 | **.49657** |

The selected intrinsic architecture is therefore width600, depth3, normalized
aggregation and a two-hidden-layer FFN.

## Confirmation split

A second connectivity split (`validation_seed=1781`) compares broad and focused
training with the same selected architecture and observation-proportional
sampling. Broad continuous training (`86416a1c`, 53 eligible heads) scores
`.47408` macro and `.39958` observation-weighted MSE. Focused CYP continuous
training (`f8c7c833`, 14 heads) improves these to `.44756` and `.37154`.

The focused objective therefore improves the common CYP evaluation suite on
both development splits. The result supports fewer, target-relevant regression
tasks; it does not support forcing a particular historical task count.

## Decision

Use 14 assay-specific CYP continuous heads, observation-proportional task
sampling, no classification or ranking objective, width600/depth3 and FFN2 as
the current project-owned pretraining recipe. Keep nine-head sharing as a
negative result. No downstream LightGBM/TabPFN evaluation, ensemble change,
candidate generation or submission was performed in this screen.

## Subsequent raw LightGBM transfer check

At the user's request, downstream evaluation was then enabled. Atomic readout
`8d226e39` and audit `555f5fe7` apply the same raw600 LightGBM protocol to the
intrinsic winner. Its macro ST-RAE is `.716760`, worse than the previous single
owned encoder `.681719` and Jeremy raw600 `.651179`.

To exclude an unlucky intrinsic selection, run `1b3e7ef9` fits the same
LightGBM readout to all six architecture-screen checkpoints and the focused
confirmation-split checkpoint. Its final aggregation failed after all seven OOF
readouts were written because a Polars Series was passed directly to
`numpy.isfinite`; the retained readouts are valid and the check is fixed in
source. Recorded audit `70e4114f` is authoritative.

The best focused checkpoint downstream is width600/depth2/FFN3 at `.706066`.
Its target ST-RAE values are `.772569/.622424/.869802/.559469` for
1A2/2C9/2D6/3A4. The best focused 3A4 value is `.538478` from the confirmation
split. No focused checkpoint beats Jeremy on any target. Depth2 improves D6
over the previous single owned encoder `.879357→.869802`, but remains worse
than Jeremy `.858244`.

Intrinsic external-task accuracy and transferable challenge representation are
therefore misaligned in this experiment. Do not replace the previous owned or
Jeremy encoders with the 14-task focused models, and do not add them to the
prediction bank. No ensemble or submission was produced.

## Downstream-guided objective sampling

Run `b5a003ef` changes the selection criterion from intrinsic assay loss to the
actual downstream raw600 LightGBM ST-RAE. It holds the single-checkpoint
width600/depth3/FFN3 encoder and endpoint-only corpus fixed, trains seven task
sampling/objective mixtures, selects independently by endpoint on challenge
fold 3, and checks only the four selected endpoint winners on untouched fold 4.
Jeremy is a frozen reference and is never used as initialization.

The fold-3 winners are source-uniform mixed for 1A2 (`.702219`), square-root
task-count mixed for 2C9 (`.536455`), square-root task-count without
classification for 2D6 (`.900233`), and observation-proportional mixed for 3A4
(`.482997`). Their corresponding fold-4 results are:

| Endpoint | Selected objective | Owned selected | Previous owned | Jeremy |
|---|---|---:|---:|---:|
| 1A2 | source-uniform mixed | .725211 | .760728 | **.705816** |
| 2C9 | sqrt-count mixed | .623230 | .607549 | **.599222** |
| 2D6 | sqrt-count, no classification | **.858356** | .875718 | .859149 |
| 3A4 | observation-proportional mixed | **.490545** | .497954 | .494655 |
| Macro | target-specific selection | .674336 | .685487 | **.664710** |

Downstream-guided task construction therefore improves the previous owned
encoder on 1A2, 2D6 and 3A4 and closes most of the macro gap to Jeremy. It does
not solve 2C9, and one universal objective is not supported: the four endpoints
select three different task-sampling designs. Continue with target-specific
source/task-family ablations. Do not change a prediction bank or create a
submission from this development screen.

## Source families and full-OOF audit

Source-family run `68f542a6` holds square-root task sampling fixed and removes
TDC, PubChem, Tox21, Polaris, ExpansionRx or PXR one at a time, plus a CYP-source
only variant. On fold 3, no-PXR appears to improve 1A2 to `.691387` and
no-Polaris appears to improve 2C9 to `.511655`. Neither result transfers to fold
4 (`.739710` and `.631860`). The no-Polaris D6 result does transfer narrowly at
`.856320`, but the effect is too small and target-specific to justify a general
source removal.

Because of this split instability, audit `7e362fe6` fits the fixed raw600
LightGBM readout across all five folds for all 14 objective/source checkpoints.
The best single owned checkpoint is square-root mixed at macro `.669477`, versus
the previous owned `.681719` and Jeremy `.651179`. Target-specific owned winners
are:

| Endpoint | Owned task design | Owned | Previous owned | Jeremy |
|---|---|---:|---:|---:|
| 1A2 | sqrt mixed, no Tox21 | .736899 | .740559 | **.728514** |
| 2C9 | sqrt mixed, no Polaris | .554558 | .597967 | **.535216** |
| 2D6 | sqrt, no classification | .870970 | .879357 | **.858244** |
| 3A4 | proportional mixed | .488025 | .508994 | **.482742** |
| Macro | target-specific owned | .662613 | .681719 | **.651179** |

Thus downstream-guided task design yields a meaningful owned-only improvement,
but Jeremy remains best on every endpoint under the complete OOF audit.

## Architecture and target-family follow-ups

Run `be329827` screens width448/depth3/FFN3, width600/depth2/FFN3 and
width600/depth3/FFN2 around each endpoint's selected task recipe. None beats the
corresponding width600/depth3/FFN3 incumbent on fold 3. Width448 is worse for all
four endpoints; the earlier focused-task depth2 D6 gain does not transfer to
the downstream-guided objective.

Run `df9fbdaa` tests an intermediate task family: retain general non-CYP ADME
tasks but keep CYP-prefixed supervision only for the target isoform. It produces
promising fold-3 values for 1A2 (`.688921`) and D6 (`.899448`) but fails to beat
the retained checkpoints on fold 4 (`.727644` and `.893357`). The 2C9 and 3A4
variants are worse on both selection and confirmation. This indicates that
cross-isoform CYP supervision is generally useful; the successful exclusions
are more specific than simply removing every non-target isoform.

Retain the four full-OOF owned winners above as research candidates. Do not
replace Jeremy, change a prediction bank or create a submission from these
screens. Further progress likely requires finer assay-level selection or new
supervision rather than more width/depth reduction.
