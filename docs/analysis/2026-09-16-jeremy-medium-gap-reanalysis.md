# Jeremy-medium gap reanalysis

## Question

Why does the externally supplied `chemprop_medium.pt` checkpoint still produce
better raw-embedding LightGBM transfer than the project-owned encoders, even
after the public source families were assembled and the owned task design was
optimized against downstream performance?

This is a diagnostic analysis. It does not train a model, change a prediction
bank or create a submission.

## What the checkpoint establishes

Direct inspection of the retained checkpoint and the installed Chemprop
implementation establishes more than the absence of task names:

- bond message passing has width 600, depth 4, no bias, no dropout and directed
  messages;
- graph aggregation is `NormAggregation(norm=100)`;
- the predictor is a 47-output `RegressionFFN` with input and hidden width 600,
  three hidden layers and no dropout;
- every output has unit task weight and a stored task-specific mean and scale;
- Chemprop's standard finite-target mask and MSE criterion jointly optimize all
  observed cells in each molecule batch; and
- the stored scheduler parameters are five warmup epochs and learning rates
  `5e-4 -> 5e-3 -> 5e-6`.

The checkpoint has no output-column names, epoch or global-step metadata, so
the complete 47-assay manifest and total optimization exposure remain unknown.
The available evidence nevertheless identifies it as a wide, standardized,
47-output numeric model trained with the standard Chemprop multitask layout.
`RegressionFFN` and MSE do not imply that every original target was continuous:
the stored transforms show five heads with Bernoulli-like means and standard
deviations, indicating that some binary targets were fitted as numeric 0/1
regression outputs.

## Recovering task blocks from the weight space

The final 47-by-600 output matrix, output-transform means/scales and the public
source tables recover part of the missing task order. This is stronger than
clustering weights alone because the saved target moments can be matched to
source-column moments.

Jeremy's checked-in `cyp_chembl.csv` is a 13,841-molecule, 12-column CYP potency
matrix whose curation follows the same assay-classification and quality-filter
logic as [Capricho](https://github.com/David-Araripe/Capricho). Eight consecutive
checkpoint heads match that table:

| Head | Checkpoint mean / scale | Recovered source column | Source mean / SD |
|---:|---:|---|---:|
| 30 | 5.419 / .856 | CYP3A4 IC50 | 5.418 / .854 |
| 31 | 5.593 / 1.172 | CYP3A4 Ki | 5.571 / 1.162 |
| 32 | 5.530 / .968 | CYP2C9 IC50 | 5.526 / .968 |
| 33 | 6.436 / 1.103 | CYP2D6 IC50 | 6.427 / 1.099 |
| 34 | 6.180 / 1.272 | CYP1A2 IC50 | 6.184 / 1.275 |
| 35 | 7.033 / 1.140 | CYP1A2 Ki | 7.063 / 1.137 |
| 36 | 5.285 / .731 | CYP2C19 IC50 | 5.285 / .726 |
| 37 | 6.576 / .936 | CYP2B6 IC50 | 6.576 / .951 |

The four omitted ChEMBL columns are the low-count Ki tasks for 2C9, 2C19, 2D6
and 2B6 (51, 136, 128 and 66 observations), consistent with an eligibility
threshold rather than arbitrary task selection.

Heads 38--41 also match the four Polaris Novartis-surrogate fields: three
pIC50 predictions centered near 4.8 and `log_kobs` centered near -1.95. Heads
42--46 have means `.439/.371/.244/.463/.459` and scales equal to the expected
Bernoulli standard deviations, strongly indicating a five-task binary block;
the five TDC Veith CYP inhibition panels are the most plausible mapping, but
their individual ordering is not asserted without an exact source snapshot
match. Output-weight cosines independently show several related task blocks,
including heads 22--26 and near-duplicate directions for heads 6/8.

Thus 12 heads are source-matched by moments, five more are identified as a
binary block, and 30 remain unnamed. Weight-space analysis cannot recover
molecule membership, missing-label patterns, sample counts or exact task names
from the checkpoint alone, but it can guide reconstruction and reject an
all-continuous interpretation.

The owned trainer is materially different. It samples one assay/objective per
optimizer step, optionally mixes continuous interval, ranking and
classification losses, and uses AdamW with a custom warmup/cosine schedule,
peak learning rate `1e-3`, final learning rate `1e-4`, weight decay `1e-4` and
dropout `0.1` in the current recipes. Consequently, matching only width, depth
or approximate task count does not reproduce the Jeremy training regime.

## Public-overlap diagnostic

Recorded run `20260916T140248000283Z-8d79b936` partitions each target's full
OOF predictions according to overlap with the currently assembled public
TDC/Polaris observations. This is only a proxy for possible upstream exposure;
it is not proof that any molecule was used to train Jeremy.

Positive deltas below mean that Jeremy is better than the target-specific owned
checkpoint.

| Endpoint | Non-overlap n | Owned - Jeremy | Public-overlap n | Owned - Jeremy |
|---|---:|---:|---:|---:|
| CYP1A2 | 1,273 | +.00821 | 99 | +.01088 |
| CYP2C9 | 1,183 | +.02984 | 72 | -.08578 |
| CYP2D6 | 1,288 | +.02224 | 169 | -.06130 |
| CYP3A4 | 2,202 | +.00726 | 98 | -.03740 |

Jeremy's advantage persists on non-overlap compounds for every endpoint. On
the small public-overlap subsets, the owned model is better for 2C9, D6 and
3A4. The only conspicuous opposite slice is same-isoform 1A2 overlap
(`n=28`, delta `+.09827`), which is too narrow and inconsistent with the other
targets to explain the global result. Bootstrap uncertainty is broad for most
slices; D6 non-overlap is the clearest difference, with a 95% bootstrap interval
of `+.00527` to `+.04182`.

The available overlap evidence therefore does not support public-data exposure
as the main explanation for Jeremy's transfer advantage.

## Representation geometry

The same run measures the three raw600 embeddings on the same 6,146 challenge
compounds.

| Embedding | Dead dims | Median feature SD | Effective rank | Participation ratio | PCs for 90% variance |
|---|---:|---:|---:|---:|---:|
| Jeremy medium | 61 | .00995 | **89.42** | **32.37** | 159 |
| Owned sqrt-mixed | 112 | .00506 | 40.29 | 13.78 | 65 |
| Previous owned | 144 | .00532 | 32.28 | 11.56 | 54 |

Linear CKA is `.7415` between Jeremy and the selected owned model, `.7239`
between Jeremy and the previous owned model, and `.8351` between the two owned
models. The owned variants therefore resemble each other more than either
resembles Jeremy and occupy a much lower-dimensional subspace.

A separate read-only activation check across 6,896 molecules reaches the same
conclusion. Jeremy has a `.215` zero-activation fraction and a median 475 active
dimensions per molecule. The selected owned encoder has `.502` and 301; the
previous owned encoder has `.472` and 320. Message-passing weight norms are also
substantially larger in Jeremy than in either owned checkpoint. The owned
models have learned away from initialization, but not to the same dense,
high-rank representation.

## Explanations ruled down

- **Depth alone:** the earlier matched depth-4 owned checkpoint scores macro
  ST-RAE `.686153`, worse than Jeremy `.651179` and the depth-3 owned checkpoint
  `.681719`. It retained the owned sampler, dropout and optimizer, so this is
  not a test of Jeremy-style training mechanics.
- **Task count alone:** reducing supervision to focused CYP tasks improves the
  held-out external task loss but worsens downstream raw LightGBM transfer.
- **Only missing public data:** Jeremy remains better in the much larger
  non-overlap slices, while owned models often win on the proxy-overlap slices.
- **A target-specific architecture tweak:** width448, depth2 and FFN2 follow-ups
  around the downstream-selected task recipes do not close the full-OOF gap.

## Leading hypothesis and next controlled test

The leading explanation is representation undertraining caused by the owned
optimization layout, rather than simply too few molecules or the wrong number
of heads. Jeremy jointly updates the encoder from every observed target cell in
a molecule batch, with standardized continuous targets, no dropout and a peak
learning rate five times larger. The owned trainer updates from one sampled
task at a time and produces a sparse, low-rank embedding even after task and
source selection improves the downstream score.

The next controlled experiment should train from scratch on the existing
challenge-disjoint owned corpus with:

- one wide target matrix and joint masked MSE over all observed cells;
- standardized numeric targets, first using the recovered eight-column
  Capricho-style ChEMBL block, four-column Polaris block and explicitly matched
  public ADME tasks; a continuous-only ablation should remain separate;
- width600, depth4, norm100, FFN3 and dropout0;
- the standard Chemprop `5e-4 -> 5e-3 -> 5e-6` schedule with five warmup epochs;
- enough complete molecule-batch epochs rather than a fixed number of
  task-local updates; and
- the unchanged target-specific raw600 LightGBM full-OOF evaluation.

A matched control should retain the same data and architecture but use the
current one-task sampler. This directly isolates the joint masked-matrix
training mechanism. It does not initialize from Jeremy, concatenate
embeddings, or require reconstructing the unknown historical 47 tasks.

## Decision

Treat the joint multitask optimization regime as the next hypothesis. Do not
attribute the gap to leakage, do not add a checkpoint to a bank, and do not
create a submission from this diagnostic alone.

## Joint-matrix experiments

Runs `20260916T143051880581Z-7c464f30`,
`20260916T143736207347Z-cc6eb40b` and
`20260916T144515536881Z-d5c8a11f` test the hypothesis from scratch on the
78,449-compound challenge-disjoint owned corpus. All variants use the same
width600/depth4/norm100/FFN3 encoder, joint masked cell-wise MSE, per-task
standardization, no dropout and the `5e-4 -> 5e-3 -> 5e-6` schedule. The
downstream comparison is unchanged raw600 LightGBM five-fold OOF.

The objective is better target-specific transfer, not reconstruction of the
unknown Jeremy training set or its 47 tasks. In particular, the experiments
include broad owned supervision and explicit ablations even when they move
away from the recovered Jeremy layout.

| Variant | Tasks | 1A2 | 2C9 | D6 | 3A4 | Macro | Effective rank |
|---|---:|---:|---:|---:|---:|---:|---:|
| all compact tasks | 72 | .74523 | .57976 | .89116 | **.48827** | .67610 | 64.01 |
| CYP-only compact tasks | 28 | .74499 | .59215 | **.85956** | .51491 | .67790 | 46.56 |
| recovered public core | 38 | .78453 | .60209 | .86775 | .52588 | .69506 | 63.64 |
| public core CYP only | 18 | .76301 | .60958 | .86340 | .52124 | .68931 | 49.27 |
| **all continuous tasks** | **59** | **.72953** | **.56196** | .86915 | .49236 | **.66325** | 47.91 |
| continuous, assay-granular | 63 | .74835 | .59073 | .88031 | .50974 | .68228 | 57.89 |
| continuous, no Polaris | 55 | .74935 | .57167 | .90238 | .53199 | .68885 | 49.74 |
| continuous, no Tox21 | 54 | .75821 | .57401 | .88476 | .49319 | .67754 | 56.89 |
| previous owned control | - | .74056 | .59797 | .87936 | .50899 | .68172 | 32.28 |
| Jeremy control | 47 | .72851 | .53522 | .85824 | .48274 | .65118 | 89.42 |

The broad continuous-only matrix is the strongest universal owned joint model.
It nearly matches Jeremy on 1A2 (`+.00101` ST-RAE) and improves the previous
owned control on every endpoint. CYP-only supervision is independently the
best new D6 representation, only `+.00131` behind Jeremy. The all-task compact
model gives the best new 3A4 result, but the earlier target-specific owned
model remains marginally better (`.48803` versus `.48827`). For 2C9, the
earlier target-specific owned model also remains better (`.55456` versus
`.56196`).

The reconstruction-like source subsets are negative. Removing all binary
tasks helps markedly even though Jeremy itself appears to contain five binary
heads. Splitting compact endpoints into assay-granular heads, removing Polaris,
or removing Tox21 all worsens the intended target as well as the macro score.
The single-task Polaris-removal hint therefore does not transfer to joint
training. Task count, named-source similarity and validation loss are not
reliable proxies for downstream representation quality.

Selecting independently across project-owned experiments gives
`.72953/.55456/.85956/.48803` for 1A2/2C9/D6/3A4, descriptive macro `.65792`.
This is better than the previous owned target-specific macro `.66261` and close
to Jeremy `.65118`, but it is selection-biased and is not an unbiased new
generalization estimate.

## Updated decision

Retain `joint_all_no_binary` as the owned 1A2 candidate and
`joint_compact_cyp` as the owned D6 candidate. Retain the earlier specialized
owned choices for 2C9 and 3A4. The experiments validate joint sparse-matrix
training as a useful mechanism but reject literal Jeremy reconstruction as the
optimization target. Do not alter a prediction bank or create a submission
from this representation screen alone.

## Target weighting, nested readout and weight-space follow-up

Runs `20260916T150348818193Z-aa2923e8` and
`20260916T151717448226Z-ea34957f` compare 2x/4x target-isoform loss weights with
matched 1x controls. Initial weights and molecule-batch order are identical
within each comparison. Every intended target worsens:

| Target | Matched 1x | 2x target tasks | 4x target tasks |
|---|---:|---:|---:|
| 1A2 | **.75686** | .76057 | .76098 |
| 2C9 | **.57413** | .57808 | .57923 |
| D6 | **.86899** | .87340 | .87120 |
| 3A4 | **.50080** | .51267 | .50283 |

The stronger weights also tend to reduce representation rank. Broad joint
supervision is useful, but directly steering its loss toward the downstream
isoform does not improve raw transfer.

Run `20260916T153455321020Z-636a0e78` then compares the four target-specific
owned raw600 candidates with Jeremy under the same honest nested LightGBM
selection. Each outer parent fold chooses one of eight readouts using mean
official ST-RAE over three inner group folds, without reading the outer labels.

| Target | Owned nested | Jeremy nested | Owned fixed | Jeremy fixed |
|---|---:|---:|---:|---:|
| 1A2 | .73401 | **.72673** | .72953 | .72851 |
| 2C9 | .55189 | **.53143** | .55456 | .53522 |
| D6 | .86029 | **.85614** | .85956 | .85824 |
| 3A4 | .49050 | **.47864** | .48803 | .48274 |

Nested tuning modestly improves owned 2C9, but worsens the other three. Jeremy
improves on all four, so readout capacity does not close the representation
gap. The earlier `20260916T152401260032Z-7c9d0237` run completed the owned half
and then stopped because its Jeremy alias pointed to OOF predictions instead
of embeddings; it is retained as a documented partial result and superseded by
the complete run.

Finally, run `20260916T160341146483Z-b76ae200` interpolates message-passing
weights at 25/50/75% between the matched-initialization 1x and 2x models. Each
interpolation remains one encoder and one raw600 embedding; it is neither an
embedding concatenation nor a prediction ensemble. No interpolation beats the
best existing owned candidate. D6 at 75% improves both matched-seed endpoints
to `.86745`, but remains worse than the earlier CYP joint model `.85956` and
Jeremy `.85824`. The other intermediate minima are 1A2 `.75594`, 2C9 `.57774`
and 3A4 `.50006`, all clearly behind their incumbents. A preceding extraction
run `20260916T160154189582Z-e73cc404` used a mismatched TF32 setting and is
explicitly discarded; the corrected run restores the training-time setting.

These follow-ups do not support a small loss-weight, readout or local
weight-space adjustment as the remaining solution. Surpassing Jeremy with one
project-owned checkpoint now likely requires a materially different
pretraining signal or training formulation, not further fine-grained tuning of
the current recipe. No bank, inference or submission change follows.

## Representation distillation

Run `20260917T051118981633Z-a9758b6c` tests Jeremy raw600 output geometry as
the missing training signal. It uses the 78,449 challenge-connectivity-disjoint
external compounds and never reads challenge pIC50, log2FC or Emax labels during
encoder training. The frozen teacher provides per-dimension standardized raw
coordinates plus within-batch cosine-similarity geometry. All students retain
the same width600/depth4/norm100 architecture and are evaluated by the fixed
target-specific raw LightGBM OOF protocol.

| Arm | Initialization/objective | 1A2 | 2C9 | D6 | 3A4 | Macro |
|---|---|---:|---:|---:|---:|---:|
| Teacher-only student | random; distillation only | **.72081** | .52595 | .86747 | .47772 | .64799 |
| Owned hybrid | owned; 59 tasks + distillation | .73533 | .57671 | .87002 | .51114 | .67330 |
| Anchored task student | Jeremy; 59 tasks + distillation | .73012 | **.51913** | .86116 | **.47498** | **.64635** |
| Frozen Jeremy | frozen control | .72851 | .53522 | **.85824** | .48274 | .65118 |

The random student beating its teacher on three endpoints rules out exact
recovery of Jeremy's 47 training tasks as necessary for transferring most of
the useful signal. It has effective rank `61.88` versus the teacher's `89.42`,
normalized challenge embedding MSE `.09065` and row cosine `.86469`, consistent
with a smoothed or compressed approximation rather than a copied embedding.

The Jeremy-initialized anchored model raises effective rank to `108.46` while
retaining cosine `.96550`. External supervision improves 2C9 and 3A4 without
fully erasing D6. Starting from the owned basin fails: its cosine reaches only
`.58886` and every endpoint is inferior to at least one control. The useful
result is therefore initialization- and representation-dependent, not evidence
that adding the current external task loss is generically beneficial.

Use the teacher-only student as the new 1A2 atomic candidate, the anchored task
student for 2C9 and 3A4, and retain frozen Jeremy for D6. This target-specific
selection is descriptive and must be handled as separate peers in later flat
bank audits; no ensemble or submission is changed here.

## PubChem-only D6 adaptation

Jeremy's public reproduction-kit README states that adapting
`chemprop_medium` on PubChem AID 1851 qHTS alone produced its best single-model
CYP2D6 result. Run `20260917T091542338422Z-b0467aeb` tests that specific signal
with the local five-head continuous matrix: 37,483 training cells overall,
including 5,038 D6 cells. No challenge labels are used.

| Arm | Teacher MSE | Cosine | 1A2 | 2C9 | D6 | 3A4 | Macro |
|---|---:|---:|---:|---:|---:|---:|---:|
| PubChem only | .98862 | .75242 | .74537 | .59880 | **.85817** | .53696 | .68482 |
| Anchor .1 | .10564 | .96060 | .73005 | .53604 | .87277 | .49419 | .65826 |
| Anchor .3 | .04787 | .98204 | **.72392** | **.51949** | .87342 | .49074 | **.65189** |
| Frozen Jeremy | 0 | 1 | .72851 | .53522 | .85824 | **.48274** | .65118 |

The unanchored arm gives the prespecified D6 target a nominal `.000075`
improvement, too small to treat as a practical gain. It also substantially
damages the other endpoints. Both anchors preserve teacher geometry but worsen
D6. The `.3` arm's 1A2/2C9 improvements overlap the stronger RL0110 candidates
and do not justify another bank member.

The reported PubChem-only direction is therefore only marginally reproduced
under this exact curation and fixed LightGBM evaluation. Jeremy's unpublished
full-research implementation may differ in PubChem processing, optimization or
downstream readout. Frozen Jeremy remains the D6 incumbent; no bank or
submission change follows.
