# Primary selection, population reversal, and D6 decision validation

The continued investigation establishes a concrete validation-population problem,
reproduces its consequence on known A4 outcomes, and tests attempted remedies.
It does not identify the realized molecular cause of the previous D6 Live losses.
The existing D6 ensemble is not adopted as a modeling base, and 49 positives are
not a constraint. No challenge payload was uploaded or replaced.

## Observed selection and retrospective connection

Run `20261002T041315546883Z-49b9e5f8` joins the supplied single-concentration
primary screen to canonical molecule identities. All 1,495 primary-covered D6
class labels lie below primary log2fc -0.5. The 2,823 development molecules at
or above -0.5 have **zero** observed D6 class labels. Two further class-labelled
molecules lack primary values and are negative. Missing class labels remain
missing throughout the investigation.

Primary-conditioned selection weights, fitted outside each held-out group and
capped at five, cannot reconstruct outcomes where there are no labels. Their
effective labelled sample is about 1,442 of 1,497. On the frozen September 15 to
September 23 AP-mass90 comparison, hypothetical prior-.07 projected UMAP MCC
delta changes from +.00708 to +.00078, still opposing the actual Live change
-.07381. Scaffold remains +.04288. Earlier successful changes retain positive
directions, but those earlier diagnostics score labelled molecules only and are
not deployment-complete comparisons. These transitions share one Live dataset;
they are not independent prospective experiments.

Run `20261002T040603473885Z-afdd90d5` compares a matched fresh neural baseline,
D6 primary auxiliary supervision, four-enzyme primary supervision, and primary
supervision plus selection weighting. Two seeds, original UMAP/scaffold groups,
five bag models and three deeper threshold folds are evaluated. The 100-epoch
baseline has mean centered MCC .09288 UMAP / .16838 scaffold. D6 auxiliary
supervision gives .07445 / .17624; cross-CYP .06189 / .17122; selection weighting
.07949 / .17254. AP improves slightly but there is no consistent decision gain.
Unused primary heads change initialization relative to earlier 16-head runs;
only this batch's matched baseline is the controlled comparator.

## Controlled A4 experiment: an actual reversal

A4 has 3,584 known classes: 1,291 are in the strong-primary region, containing
605 positives; 2,293 are weak-primary or primary-unavailable, containing 159
positives. Unlike D6, these provide known outcomes for a masking experiment.

Runs `20261002T042725116238Z-fe96ab3c` (A4) and
`20261002T042715697470Z-2f48e876` (D6) compare matched 20-epoch Morgan networks:
all known classes, strong-primary classes only, and the same selected-class fit
with its predicted primary-activity probability multiplied into its class score.
Activity means measured primary log2fc below -0.5; it is not the TDI definition.
Class losses and threshold selection both exclude hidden weak/missing-primary
classes. Outer evaluation restores their actual labels. Molecular groups and all
assays are excluded from every held-out fitting boundary. Primary activity targets
are measured-screen targets; no unlabelled TDI molecule becomes a negative.

Mean results over two seeds for all-class versus selected-class training:

| Evaluation population | Split | AP all -> selected | MCC all -> selected |
| --- | --- | --- | --- |
| Strong primary only | UMAP | .47839 -> .53528 | .03060 -> .07544 |
| All known A4 classes | UMAP | .29017 -> .22018 | .16330 -> .00568 |
| Strong primary only | scaffold | .49994 -> .55550 | .04181 -> .14617 |
| All known A4 classes | scaffold | .34812 -> .23697 | .21445 -> .05455 |

Selected-minus-all MCC on the full known population is -.15762
[-.23429,-.07614] UMAP and -.15990 [-.19568,-.12248] scaffold, conditional group
bootstrap intervals. On selected scaffold classes alone the direction reverses:
+.10437 [.04291,.16856]. The strong-only gain is not a gain on the broad population.

The reversal survives matched full-pool rank fractions; it is not solely caused by
different fitted positive counts. At a 7% full-pool rank fraction, mean A4 MCC is
.06358 -> .00086 UMAP and .13147 -> .02247 scaffold. This is a controlled analogy
to fixed-count ranking failures, not an assertion of the unknown D6 blind labels.

The primary-activity product does not rescue the selected learner: full-known A4
MCC .00710 UMAP / .04305 scaffold, and D6 decisions deteriorate as well. A primary
cutoff cannot be assumed necessary for TDI: the masked A4 region contains 159
actual positives. D6 weak-region metrics cannot be estimated; the two negative
primary-unavailable labels are not representative of its 2,823 weak molecules.

## Assay-compatible gate and unconstrained thresholds

Runs `20261002T043025957609Z-a3d07bb8` (D6) and
`20261002T043105536476Z-8e0b5aca` (A4) test a separate condition compatible with
the observed class definition: TDI-arm pIC50 > 4+log10(2). Every known positive
passes it. A three-head Morgan network learns class, measured potency activity,
and primary regression. Conditional class training restricts class loss to
potency-active molecules; threshold fitting still uses all known training classes.
The gate multiplies conditional class and predicted potency-activity probabilities.
It is an estimated probability decomposition, not a calibrated joint posterior.

A4 improves under this decomposition; D6 is seed-dependent. D6 gated MCC is
.10980/.08079 UMAP and .13156/.14680 scaffold. Conditional training alone is
.10722/.05073 and .13527/.13842. These do not establish a reliable new D6 winner.
Neither physical consistency nor adding primary assays resolves missing coverage.

Run `20261002T043618387974Z-bc88d8eb` tests genuinely train-selected thresholds
at illustrative positive priors .03, .07 and .12, using identical fits across each
policy. Counts are unconstrained. At matching held-out prior .07, prior-trained
cuts versus ordinary centered cuts change projected MCC by -.00032
[-.02106,.01759] UMAP and -.03884 [-.06122,-.01518] scaffold for the fresh class
network. Potency-gate changes are -.00764 / -.03909. Changing the assumed prior
does not establish a decision gain. These scenarios assume label shift and are
not identified blind accuracy; the selection gap violates the needed coverage.

Standalone activity/potency controls used Torch's highest matrix precision;
the mixed fresh-network/prior batch used high precision after the fresh fit set
it. Comparisons are within batches, not exact cross-batch reproductions. Source
and configs now set these choices explicitly per model to remove model-order
dependence from subsequent runs.

Run `20261002T044807757632Z-779faca0` evaluates the remaining observed seed
instability by equal averaging of the two saved seed margins. It first exactly
reconstructs every original call from the saved mean bag cut; no outer label
enters a fit, threshold, or ensemble weight. The resulting ten-model fresh
network has MCC .12572 UMAP / .16871 scaffold, versus individual-seed means
.11844 / .16513. Paired deltas +.00727 [-.01003,.02387] and +.00358
[-.01559,.02195] are small. Potency-gate averaging improves UMAP but worsens
scaffold and remains below this fresh class learner in both splits. Averaging
addresses random fit variation without establishing a material decision gain or
recovering the absent weak-region class outcomes. The proposal follows observing
seed instability; conditional intervals do not cover that protocol selection.

## Completed decision and records

The experiments support a mechanism that makes selected-population AP and MCC
misleading, with an observed A4 reversal and a D6 coverage gap. They do not show
that this caused the precise historical D6 swaps or that an existing D6 model is
valid. Neither selection weighting, the two gates, nor prior-based thresholds
establishes a new D6 gain worth an improvement-seeking twelve-hour upload.
A Live score would evaluate a combined pipeline without revealing the missing
weak-region labels. No new submission is recommended from this completed batch;
this conclusion does not reinstate the incumbent or a 49-call budget.

Reviews `20261002T043637506788Z-2f9a88ee` and its extended matched-rate companion
`20261002T043821637258Z-301d0e01` preserve population metrics, paired group
intervals, seed agreement, and full-pool comparisons. Prior review
`20261002T044003704380Z-5a9fd225` preserves matching-prior metrics and intervals.
Intervals condition on saved fits and omit training and protocol-selection
uncertainty; original development labels have already been inspected.

Initial A4 run `20261002T042640901914Z-28092609` was stopped after detecting
incomplete masking in the single-refit comparator. Its partial outputs and reason
are preserved and excluded from comparisons. The corrected complete run is used.
All eleven completed records are promoted to `artifacts/tracked/cyp-d6-zero-start-*`
with source-run IDs, input hashes and configs. Recorded experiments, fitting-boundary
assertions, `pixi run check`, and 73 passing / 8 skipped infrastructure tests
provide the implementation checks. The earlier 77-call prepared payload remains
immutable and unsent.
