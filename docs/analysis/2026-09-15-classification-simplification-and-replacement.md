# Classification ensemble simplification and peer replacement

## Decision

Keep the currently scored D6/A4 boolean vector unchanged. The development
system can be made easier to reason about, but the evidence supports different
actions by endpoint:

- D6: retain the current hierarchy as the submission control. A common OOF
  prediction bank is now the development interface, but neither unrestricted
  learned stacking nor anchor-preserving stacking beats the control's
  fixed-rate MCC. The equal mean of all 24 experts improves AUROC/AP while
  tying MCC and is the cleanest non-heuristic research alternative.
- A4: reduce the development Caruana pool from seven top-level candidates to
  five by removing `lgbm::morgan_balanced` and
  `pfn::external_class_pca256`. The five-member pool preserves cross-fit
  performance; the three-member pool does not.
- Frozen CheMeleon embeddings: use raw parent embeddings as the primary A4
  readout and PCA only as a compression control. PCA is not a generally
  beneficial preprocessing step here.
- SMARTCyp, ATTNSOM, CYP-MAP, scored BioTransformer products, raw
  parent/product embeddings and MACE remain representation-bank experiments.
  Adding all of them to Caruana increases selection freedom without improving
  held-out performance.

This analysis creates no candidate, submission, or upload.

## What the current scored model actually contains

The two endpoints do not share an ensemble design.

### CYP2D6

The deployed test score is hierarchical:

1. `old11`: the equal probability mean of six LightGBM models, three TabPFN
   models, and two graph models.
2. `live = rank(2 * old11 + Morgan-logistic) / 3`, where the notation means a
   2:1 rank blend rather than a probability-scale average.
3. `deployment = rank(2 * live + ComplementNB) / 3`.
4. Select the top 49 of 750 compounds.

The honest OOF estimate is not obtained by applying ComplementNB in every fold.
The historical second-stage screen chose `live+RBF-SVM` on held folds 0 and 1
and `live+ComplementNB` on held folds 2--4; full-OOF deployment chose
ComplementNB. This distinction explains why the correct current cross-fit
target-prior MCC is `.09988`, not the `.11011` obtained by retrospectively
using the full-fit ComplementNB choice in every fold.

### CYP3A4

The current score is an 80/20 rank blend of:

- a six-effective-member AUROC-Caruana prediction: historical
  CheMeleon-PCA256 + nested predicted-log2fc4 + Mordred TabPFN (`.48`), Jeremy
  PCA256 TabICL (`.20`), Jeremy PCA256 TabPFN (`.16`), external-class+Mordred
  balanced LightGBM (`.10`), external-class PCA256 TabPFN (`.04`), and Morgan
  balanced LightGBM (`.02`);
- Morgan + five classical PXR pEC50 predictions + four frozen Admetica CYP
  inhibition probabilities (`.20`).

The top 249 of 750 compounds are called positive. SMARTCyp, ATTNSOM, CYP-MAP,
BioTransformer, MACE and strict PXR-Chemprop are not in this scored model.

## Why earlier decisions looked inconsistent

Three different questions had been mixed together:

1. **Standalone representation value:** does a feature family improve a weak
   Morgan readout on UMAP and scaffold splits?
2. **Fixed incremental value:** does a prespecified 10% or 25% rank addition
   improve the current live score?
3. **Peer/replacement value:** if old and new candidates are placed in the same
   selection pool, does fold-excluded Caruana retain the new family or use it to
   replace an incumbent?

Most SMARTCyp/ATTNSOM/CYP-MAP and BioTransformer work had answered (1) and (2),
not (3). Strict PXR-Chemprop was the exception: it had already been generated as
a fully nested peer and was rejected by outer-fold fixed-rate evaluation
(`2e6c1f38`, `9a7f115d`). The new peer audit closes this gap for the other saved
families with second-level cross-fitting over their existing OOF predictions.
It is less strict than refitting every base model inside another nested loop,
but held-fold labels never influence member selection.

## Raw versus PCA CheMeleon embeddings

Run `20260915T073340818290Z-d781ddd2` compares Morgan plus raw 2048-dimensional
CheMeleon features with fold-local PCA32/128/256. `parent` is the parent
embedding, `delta` is the endpoint-specific mean product-minus-parent vector,
and `parent+delta` concatenates the two.

### A4 target-prior MCC

| View | UMAP raw | UMAP best PCA | Scaffold raw | Scaffold best PCA |
|---|---:|---:|---:|---:|
| Morgan only | .18788 | -- | .31435 | -- |
| Parent | **.31745** | .27548 (128) | **.38241** | .36698 (128) |
| Product delta | .21456 | **.22552 (128)** | .27237 | **.30814 (128)** |
| Parent + delta | **.28482** | .28482 (128) | **.35307** | .34688 (128) |

Raw parent features are decisively preferable: AUROC/AP are also
`.7323/.4336` on UMAP and `.7837/.4769` on scaffold, versus
`.7105/.3948` and `.7746/.4624` for PCA128. The metabolite delta is much weaker
than the parent and sometimes benefits from PCA128, consistent with noisy
high-dimensional differences. Concatenating the delta does not beat the raw
parent.

### D6 target-prior MCC

Morgan alone scores `.04284/.10153` on UMAP/scaffold. No embedding view improves
both splits. Parent PCA256 only ties Morgan on UMAP and falls to `.01472` on
scaffold; raw parent scores `.02959/.02743`. The best isolated scaffold delta
result (`.07710`, PCA256) still trails Morgan. D6 therefore supplies no reason
to adopt either raw or PCA embeddings.

PCA had a practical rationale: reduce 2048 features for small-data readouts,
regularize noisy product deltas, and make TabPFN-compatible inputs. It is not an
information-preserving default. For tree readouts that accept the raw width,
the matched result says to test raw first.

## Compression

Prediction-fidelity run `20260915T074431941749Z-dd47364d` fits nonnegative
weights to reproduce the current rank. It is a diagnostic, not independent
validation.

For A4, a three-member approximation (historical, Jeremy TabICL, PXR/Admetica)
has blind-test Spearman `.99436` to current and ten top-249 swaps, but OOF MCC
drops `.42000→.40698`. A five-member approximation adds Jeremy TabPFN and
external-class+Mordred, reaches Spearman `.99761`, seven swaps, and OOF MCC
`.41421`. Morgan receives effectively zero weight and external-class TabPFN
only `.0022` in the unrestricted fit.

The authoritative cross-fit pool comparison is run
`20260915T080727072715Z-11c62728`:

| A4 pool | AUROC-Caruana AUROC | AP | target-prior MCC | Group-bootstrap MCC delta vs 7 |
|---|---:|---:|---:|---:|
| Current seven-candidate pool | .81149 | .53157 | .45927 | reference |
| Five candidates | **.81159** | .53078 | .45467 | -.00455 `[-.01134,.00191]` |
| Three candidates | .80976 | .52526 | .43780 | -.02127 `[-.04236,-.00204]` |

With the sparse-MCC selector, the five-member pool is exactly tied with the full
pool at `.44394` target-prior MCC and has a near-zero bootstrap delta; the
three-member pool again loses. Thus five is the defensible development-pool
compression point. It removes two weak selection degrees of freedom, not the
underlying historical-stack submodels. The existing scored vector should not be
changed merely to enact this code-level simplification.

D6 is only three candidates at its final selection layer, although its `old11`
dependency is broad. Expanding and simplex-distilling all 13 underlying ranks
failed to reproduce the current score/rank reliably. This means the legacy
prediction cannot be losslessly simplified; it does not mean the same hierarchy
should remain the development interface.

## Peer-pool and explicit replacement results

The final audit uses one prespecified representative for SMARTCyp, ATTNSOM,
CYP-MAP, all-mechanism fusion, product-activity summaries, raw parent,
raw parent+product delta, and MACE. For D6 these are compared at the correct
second-stage level as `2 * live + candidate`; for A4 they enter beside the seven
top-level candidates.

| Endpoint | Selector | Existing pool MCC | Expanded pool MCC | AUROC change | MCC 95% group interval |
|---|---|---:|---:|---:|---:|
| D6 | AUROC Caruana | .11011 | .10496 | -.00654 | `[-.02532,.01161]` |
| D6 | sparse target-MCC | .11011 | .11011 | +.00258 | `[-.01359,.01096]` |
| A4 | AUROC Caruana | .45927 | .45161 | -.00073 | `[-.01966,.00296]` |
| A4 | sparse target-MCC | .44394 | .43166 | -.00251 | `[-.02852,.00186]` |

The D6 AUROC selector gives nonzero full-fit weights to ComplementNB, SMARTCyp,
all-mechanism, product activity and MACE, but choices vary substantially by
held fold and the aggregate result is worse. The A4 AUROC selector repeatedly
selects SMARTCyp at small weight, yet also performs worse than the existing
pool. Selection frequency is therefore not evidence of incremental value.

The best-looking explicit D6 swap replaces `live+ComplementNB` with raw parent
embedding and scores `.11531`, but its AUROC falls to `.57157`; the grouped MCC
delta is only `.00522 [-.01774,.02920]` with 66.3% improvement probability.
This is not robust after a large swap screen. On A4, many nominal swaps exactly
equal the base result because the removed Morgan member and the inserted family
both receive zero weight. The largest nontrivial sparse replacement
(Jeremy TabICL to SMARTCyp) reaches `.45007`, below both the full sparse pool and
the current fixed score, with an interval crossing zero.

## Genuinely flat D6 follow-up

The hierarchy-level analysis above still wrapped each new D6 representation as
`2 * live + candidate`. The user correctly noted that this preserves the same
structural problem under different names. Corrected run
`20260915T083537675615Z-6631a8ca` therefore removes both the `old11` mean and
all `live/live2` wrappers. It places the eleven historical atomic predictions,
Morgan logistic, RBF-SVM, ComplementNB, and eight new representation models in
one Caruana pool. All joins explicitly preserve the left-side molecule order
before conversion to NumPy arrays.

| D6 design | AUROC | AP | target-prior MCC |
|---|---:|---:|---:|
| Current honest legacy OOF | **.57933** | **.29059** | **.09988** |
| Flat legacy 14, sparse-MCC selection | .54369 | .26062 | .07049 |
| Flat expanded 22, sparse-MCC selection | .56271 | .27398 | .07049 |

The expanded flat pool improves AUROC over its flat predecessor by `.01893`,
but the grouped interval crosses zero `[-.00120,.04076]`; MCC is unchanged
within noise (`-.00072 [-.04173,.03691]`). Both flat designs are substantially
worse than the current hierarchy. The flat expanded full fit selects external
class+Mordred, Morgan logistic, ComplementNB, raw parent and MACE, but this
weight vector does not rescue fold-excluded performance.

The flat interface is mechanically extensible, but it is not a supported
successor architecture. The earlier run `94f83b39` that suggested otherwise is
invalid: a Polars join did not preserve molecule row order before array
conversion. Its numerical claims are withdrawn.

## D6 `old11` basic-feature ablation

Run `20260915T083537675615Z-440b1372` holds the downstream hierarchy fixed and
removes the plain Morgan, Mordred and FP+Mordred LightGBM members from the
eleven-model probability mean.

| Variant | old11 members | AUROC | AP | target-prior MCC | test swaps each way |
|---|---:|---:|---:|---:|---:|
| Full `old11` | 11 | .57933 | .29059 | **.09988** | 0 |
| Drop Morgan | 10 | .57783 | .29321 | .09485 | 2 |
| Drop Mordred | 10 | .58045 | .29183 | .08495 | 2 |
| Drop FP+Mordred | 10 | .58054 | .29230 | .08495 | 2 |
| Drop all three | 8 | **.58070** | **.29534** | .08495 | 7 |
| Drop all six LightGBM members | 5 | .58319 | .28747 | .07049 | 10 |

The three basic models are not necessary for rank discrimination: removing all
three slightly improves AUROC and AP. They do affect the fixed historical
operating recipe. The grouped MCC delta for dropping all three without changing
anything downstream is `-.01489 [-.04016,.00810]`. This is a reproduction
diagnostic, not the correct pruning decision, because a simplified model should
be allowed to choose a new downstream recipe.

## Performance-first D6 simplification

Run `20260915T105041333820Z-ba8c651c` removes test-vector fidelity from the
objective. For every pruned core it evaluates fixed simple recipes and
held-fold-excluded recipe/Caruana selection. The cleanest successful design is:

1. remove plain Morgan, Mordred and FP+Mordred LightGBM from `old11`;
2. average the remaining eight models as one core block;
3. retain the Morgan logistic block; and
4. replace the fold-specific RBF/Complement choice with their fixed average.

Equivalently, the final score is one fixed four-block blend: core `4/9`, Morgan
logistic `2/9`, RBF `1/6`, and ComplementNB `1/6`. There is no fold-dependent
deployment rule.

| Design | Core members | AUROC | AP | target-prior MCC |
|---|---:|---:|---:|---:|
| Current hierarchy | 11 | .57933 | .29059 | .09988 |
| Full core + fixed tail average | 11 | .58071 | .29445 | **.11531** |
| Drop Mordred + fixed tail average | 10 | .58179 | .29579 | .11011 |
| Drop FP+Mordred + fixed tail average | 10 | .58189 | .29548 | .11011 |
| Drop all three basic models + fixed tail average | **8** | **.58216** | **.29759** | .10496 |

The reduced design improves all three point metrics over the current hierarchy.
Its grouped MCC delta is `.00522 [-.02215,.03471]`, so the improvement is not
statistically secure. This was only a prespecified ablation of the three models
named in the discussion, however, and is not a valid claim that those are the
best members to remove.

## Exhaustive all-core simplification

Run `20260915T111425299454Z-bea02fa4` corrects that remaining scope error. It
searches every one of the 2,047 non-empty subsets of all eleven core members,
combined with five downstream recipes. Subset and recipe choice use only the
other four folds; results below are the held-fold predictions. Component count
includes the core members and any logistic/RBF/Complement blocks used.

| Components | Selection objective | AUROC | AP | target-prior MCC |
|---:|---|---:|---:|---:|
| Current 14-component hierarchy | historical | .57933 | .29059 | .09988 |
| 5 | AUROC | .58785 | .29622 | .08008 |
| 6 | MCC | .57550 | .28160 | .11011 |
| 9 | MCC | .57918 | .29392 | .10496 |
| **10** | **MCC** | **.58059** | **.29623** | **.11011** |
| 13 | AUROC | .58649 | .29837 | .11011 |
| 14 | fixed tail average | .58071 | .29445 | .11531 |

Ten components are the smallest strict point-metric Pareto replacement: nine
components miss the current AUROC by `.00016`, while ten improve AUROC, AP and
MCC. Its MCC delta is `.01093 [-.00987,.03486]` with 85.3% positive-bootstrap
probability.

The full-data ten-component choice uses a seven-member core—plain Morgan,
Mordred, FP+Mordred, external-class, external-mixed, external-mixed TabPFN and
the TDI-only graph—followed by the fixed logistic/RBF/Complement tail average.
Four of five outer folds choose this same seven-member core and recipe; the
remaining fold chooses an eight-member core with the Complement tail. Thus the
original three basic models are retained by the comprehensive search. The four
legacy members excluded from the full-data choice are external-class+Mordred,
external-class TabPFN, CheMeleon TabPFN and the auxiliary graph.

Correlation audit `20260915T112254509486Z-0ed6e4b8` shows genuine redundancy
inside the retained Morgan block. Pairwise OOF Spearman correlation among plain
Morgan LGBM, logistic, RBF and ComplementNB ranges from `.618` to `.880`; the
largest is logistic versus RBF. Their fixed-rate call Jaccard is lower
(`.190–.516`), so they still disagree near the decision boundary. Despite its
feature name, FP+Mordred LGBM correlates only `.256–.285` with the four pure
Morgan models and behaves as a distinct predictor. Descriptive leave-one-out
removal of RBF improves AUROC/AP but lowers MCC `.12058→.11011`. Thus the
Morgan block is overrepresented globally, but the present OOF evidence does not
identify a cost-free removal at the top-call operating point.

Raw external-mixed TabPFN had not been tested previously. Run
`20260915T112944344816Z-69abecf5` fits the same TabPFN2.6 recipe on scaled raw
600-dimensional embeddings without PCA. Relative to PCA256, raw improves D6
AUROC `.56065→.56736` and AP `.25882→.26835`, while target-prior fixed-rate MCC
changes `.05646→.05188`; prediction Spearman is only `.484`. Expanded exhaustive
run `20260915T113122191419Z-c156ae61` selects raw600 in every size-ten outer
fit, confirming useful signal, but the enlarged selection space yields
`.58297/.29473/.09988` and does not supersede the legacy-only ten-component
cross-fit result. Retain raw600 in the representation bank.

Controlled run `20260915T113718060068Z-b3a7120e` then holds the original
ten-component selections fixed. Replacing PCA256 with raw600 changes cross-fit
AUROC/AP/MCC from `.58059/.29623/.11011` to `.58051/.29635/.10496`.
Averaging raw600 and PCA256 as one conceptual embedding block gives
`.58076/.29686/.11011`: a small rank improvement with identical MCC. This
requires two TabPFN inference models. By user decision, use the averaged block
as the first eleven-model candidate and retain PCA256 alone as the strict
ten-model lower-cost fallback.

## Extensible learned stacking

Run `73198c15` replaces manual placement and fixed blend weights with one common
OOF prediction-bank interface. The current fixed-ratio prediction remains an
eligible expert and the control; adding a model means appending one keyed OOF
column. The bank contains 24 D6 experts and 16 A4 experts, including the legacy
atomic predictions, SMARTCyp, ATTNSOM, CYP-MAP, mechanism fusion,
BioTransformer product activity, parent/product embeddings, MACE and D6 raw600
TabPFN.

Every meta-model is evaluated by outer-fold cross-fitting. Ridge regularization
and convex-simplex penalties are selected only inside the outer training folds.
The comparison also includes a pairwise-ranking simplex and a permutation-
invariant expert-set network whose shared gate accepts a variable number of
experts. Thus the comparison evaluates an addition procedure, not a
retrospectively chosen blend.

| Endpoint | Method | AUROC | AP | target-prior MCC |
|---|---|---:|---:|---:|
| D6 | Current fixed ratio | .57926 | .29124 | **.09988** |
| D6 | Equal mean, all 24 | .58296 | .30282 | **.09988** |
| D6 | Simplex log loss | **.58510** | **.30517** | .08495 |
| D6 | Ridge | .55634 | .27280 | .05646 |
| D6 | Expert set | .54473 | .26297 | .06576 |
| A4 | Current fixed ratio | .80713 | .51273 | .45467 |
| A4 | Ridge | **.80924** | **.53012** | **.45774** |
| A4 | Pairwise simplex | .80839 | .52899 | .44854 |
| A4 | Expert set | .79996 | .51823 | .41014 |

The A4 ridge point estimate improves all three ranking/operating metrics, but
its grouped MCC delta is only `.00340 [-.02034,.02679]`; retain it as a
development candidate rather than a proven replacement. D6 simplex methods
improve global ranking but move the top-call boundary in the wrong direction.
The equal mean ties the control MCC with a near-zero grouped difference
`.00024 [-.02877,.03118]`; it is structurally clean but not evidence of a fixed-
rate improvement. The expert-set network is rejected in its first formulation.

Safe follow-up `07a0d940` penalizes a pairwise simplex toward the current expert
and separately learns only the interpolation between current and the mean of
all other experts. The D6 anchored simplex keeps an average 73% weight on the
current expert but reaches only `.58088/.28519/.08987`; the learned mean
interpolation reaches `.58069/.29379/.07526`. A4 anchor interpolation ties the
current MCC and does not improve AUROC. Anchor preservation therefore does not
rescue D6, and no learned stacker supersedes its fixed-ratio submission control.

## Simplified working policy

1. Report standalone dual-split value first; use raw frozen embeddings whenever
   the readout can handle them and keep PCA as a matched control.
2. Promote only one representative per correlated feature family into ensemble
   testing.
3. Compare `existing pool` and `existing + representatives` using fold-excluded
   selection. Report a replacement only when its bootstrap contrast is stable,
   not because it appears in a full-OOF weight vector.
4. Register new OOF predictions in the common prediction bank and rerun the
   prespecified outer-cross-fit stackers. Do not manually decide which legacy
   block receives a new model.
5. Keep the current D6 fixed ratio as the submission control and the 24-expert
   mean as a clean research alternative. Keep A4 ridge as a development
   candidate, subject to full-data inference and prospective validation.

## Run disposition

- `d781ddd2`: retained raw/PCA evidence.
- `dd47364d`: retained prediction-fidelity diagnostic.
- `efe7b479`: discarded; first D6 peer screen used the wrong hierarchy level.
- `c3d63c8c`: failed; probability-disabled SVC was called with
  `predict_proba`.
- `42ab66b9` and `abf0e2e7`: corrected intermediate audits, superseded by the
  final run after adding replacement bootstrap and compact A4 pools.
- `11c62728`: authoritative A4 compression and hierarchy-level replacement
  audit; D6 architecture conclusion superseded by the flat follow-up.
- `94f83b39`: discarded; non-preserved join order misaligned prediction arrays.
- `73156203`: failed reconstruction that exposed the join-order issue.
- `6631a8ca`: authoritative corrected hierarchy, replacement and flat audit.
- `440b1372`: authoritative fixed-policy D6 `old11` LightGBM ablation.
- `ba8c651c`: limited performance-first ablation, superseded for member choice.
- `bea02fa4`: authoritative exhaustive all-core D6 simplification audit.
- `a364728b`: interrupted during excessive small-matrix thread overhead; no
  result was produced.
- `e1b8e294`: interrupted after identifying repeated pairwise-matrix
  construction; retained as a failed implementation run.
- `73198c15`: authoritative extensible-stacking comparison.
- `07a0d940`: anchor-preserving follow-up; discarded as a replacement and
  retained as negative evidence.
