# Focused regression program: 1A2 and 2D6

The user prioritizes improving CYP1A2 and CYP2D6 after the latest LB result,
without forbidding other CYP or classification auxiliary tasks. The current
submitted five-component OOF ensemble remains the comparison reference:
1A2 ST-RAE0.683335; 2D6 ST-RAE0.834741. The official baseline is a comparison,
not the success target. No new numeric stopping goal or submission was requested.

## Prospective experiments

1. Saved prediction diagnosis and exploratory calibration across folds. Completed
   in `20260912T082119640472Z-88efb7f8`; simple calibration has little effect, while
   component residual correlation is high. Details: `2026-09-12-focus-diagnosis.md`.
2. Target-focused multitask interval training from the same upstream medium
   encoder. pIC50 weights1A2/2C9/2D6/3A4=2/.5/2/.5, log2fc weights.25 versus1.
   Checkpoint selection uses inner validation ST-RAE on1A2/2D6 only. Other labels
   remain auxiliary supervision. Evaluate direct predictions and LGBM readout.
3. Pure single-target interval adaptation for1A2 and2D6 separately. Same parent
   outer folds, grouped inner checkpoint selection, max160epochs/patience20.
   Evaluate direct prediction and embedding+log2fc LGBM to distinguish the head
   from the learned representation.
4. Two RBF-SVR readouts (frozen upstream600 or47scores, each+log2fc4) and one
   combined LGBM view (600+47+Morgan2048+log2fc4). Each uses16 Rustuna trials on
   a grouped inner holdout, then refits on the outer training labels. All five
   folds and both focus endpoints;480 fits total. Frozen features have no
   challenge-label fitting. Auxiliary features retain outer-safe cross-fitting,
   but no additional inner encoder fit; inner selection optimism is explicit.
   Run `20260912T082257132357Z-bf656b4e`.
5. Alternative primary-screen featurizer: frozen upstream features → LGBM on
   actual measured log2fc, with the same ten pair-excluded training sets. This
   tests whether the current primary predictor limits downstream information.

All training comparisons preserve the retained parent assignments and exclude
all outer compound labels from challenge neural fitting. New prediction
combinations use prespecified25%/50% additions as diagnostics. Do not silently
replace the submitted model based on a small development OOF fluctuation.
All runs, including failed attempts, remain recorded. Classification is allowed
where useful; it is not inserted solely to satisfy a model naming convention.

## Evidence-triggered auxiliary replacement follow-up

The new LGBM primary featurizer improves actual log2fc prediction on all four
heads, and an identical downstreamLGBM comparison improves all four pIC50
endpoints. Therefore evaluate replacement of the auxiliary feature source in
three existing final-ensemble components: upstream47scoreTabPFN,
CheMeleonPCA256TabICL, and original fine-tuned-mediumLGBM with the previously
selected per-fold parameters (no new tuning). Keep other features and model
parameters fixed. Compare each replacement separately and all three together in
the original equal five-component ensemble. This is a development comparison,
not a claim of independent untouched validation or a new submission decision.

## Completed result and retained candidates

All launched comparisons are finished. One optional matching-logfc single-target
arm was interrupted during coordination; its partial record remains and is not
negative model evidence. Other failed audit/early-input attempts remain recorded
in the component reports. No production model pointer or submission changed.

The useful direction was improving the actual primary-screen featurizer, not
simply increasing the focus endpoints' loss weights. Replacing its model with
LGBM on external pretrained600emb+47scores improves measured-logfc prediction
on all four CYPs. Matched downstream LGBM and TabPFN controls both improve;
the original tunedFT component also improves with the exact same20 selected
parameter sets. The auxiliary anchor inclusion policy differs from the old
Chemprop auxiliary, as detailed in the primary report.

The primary replacement comparison `20260912T083535835184Z-2971db50` preserves
all6,384 OOF rows and reconstructs the original ensemble within1e-12. Replacing
its three auxiliary-using components gives macro0.587787 (neighbors0.676836)
versus original0.596414 (neighbors0.687208). Adding the new frozen-medium+
primaryLGBM readout as an equal sixth component gives macro0.587892 and slightly
better focus-target scores,1A2 .670809 and2D6 .826527.

A separate independent readout screen found small complementary2D6 signal from
600emb+47scores+Morgan2048+oldnestedlog2fc → tunedLGBM. Retaining this at previously
screened25%/50% weights only for2D6 produces two OOF-only candidates:
`20260912T083627802028Z-ae562a59`. Other three endpoint predictions are checked
bitwise unchanged from the updated six-model blend.

| Endpoint | Submitted recipe OOF | Updated six +25% 2D6 addition | Updated six +50% 2D6 addition |
|---|---:|---:|---:|
| 1A2 | 0.683335 | 0.670809 | 0.670809 |
| 2D6 | 0.834741 | 0.822792 | 0.822681 |
| 2C9 | 0.463781 | 0.452634 | 0.452634 |
| 3A4 | 0.403800 | 0.401599 | 0.401599 |
| Macro | 0.596414 | 0.586958 | 0.586931 |

The two2D6 addition weights perform almost identically; this is not evidence to
fine-tune their continuous weight. Both are retained, with no automatic adoption.
These are development OOF improvements after comparisons, not new LB scores;
the last submitted LB remains0.721167. The gains are modest relative to the gap
that motivated this work. Classification auxiliaries were allowed but not forced;
this round does not rule them out or claim a classification improvement.

## Reproduction and evidence

- [Saved-prediction diagnosis](2026-09-12-focus-diagnosis.md).
- [Focused multitask controls](2026-09-12-focused-multitask.md).
- [Single-target adaptation](2026-09-12-focus-single-target.md).
- [Kernel/tree readouts and conditional group bootstrap](2026-09-12-focus-readouts.md).
- [Primary featurizer and paired readouts](2026-09-12-focus-primary.md).
- [Fixed-parameter FT component replacement](2026-09-12-focus-ft-primary.md).

Candidates under `artifacts/runs/regression/20260912T083627802028Z-ae562a59/outputs/`
are `focus_2d6_combined_25pct.parquet` and `focus_2d6_combined_50pct.parquet`.
The replacement run2971db50 contains each intermediate composition and its exact
input/component configuration. All encoder checkpoints, nested primary models,
readout models, selections and memberships remain saved in their source runs.
Checks pass; tests37passed/7skipped. No test inference or new upload was performed.
