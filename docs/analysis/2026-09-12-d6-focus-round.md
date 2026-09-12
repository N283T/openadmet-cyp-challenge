# CYP2D6-focused follow-up to the community working notes

The user authorized further focused experiments while away. The aim is to improve
CYP2D6 regression, retaining the existing parent folds and other CYPs when useful.
No new leaderboard submission is performed by this experiment round.

## What the current community entry does

The [current working notes](https://supercowpowers.github.io/workbench/blogs/cyp_challenge/)
add two CYP2D6 specialists to a broader multitask ensemble: one sees all available
2D6 readouts and one sees only direct 2D6 pIC50. Public assays remain separate
heads. Their maximum-response signal keeps inactive compounds in training.
The displayed entry also includes per-target placement tuned using leaderboard
feedback; its score is not an uncalibrated model-only comparison.

Source snapshot `7a70c5e3` preserves the current blog and linked scripts. The
published union configuration has evolved since the entry: its script now has
26 heads, versus 18 in the reported ensemble member. The public task weight is
per observed label, not the fraction of total training loss.

Our earlier curation already retained maximum response. The new comparisons
therefore test how to use that signal and whether specialist training helps,
rather than describing the same public data as newly acquired.

## Experiments

1. Matched full-CYP versus 2D6-auxiliary versus direct-only graph adaptation,
   all using a 2D6 inner checkpoint criterion. Fresh initialization and an
   appropriate scratch learning rate are additional controls.
2. External maximum-response, observed high-dose response and inactive prediction
   models as frozen features. See [efficacy comparison](2026-09-12-d6-efficacy.md).
3. Actual joint graph learning with separate challenge and public response heads,
   including no-public, 0.05 and 0.30 per-label public weights.
4. Inner-cross-fitted interval-width-conditioned residual decisions. See
   [width comparison](2026-09-12-d6-interval-width.md).
5. Sparse count-fingerprint kernel SVR and matched point/interval kernel
   regression. See [kernel comparison](2026-09-12-d6-kernels.md).

Completed OOFs were evaluated against the frozen previous global pool, with
both full-OOF fitted weights and held-fold weight selection reported separately.
Underlying OOF models are reused; weight cross-fitting checks selector stability.
Unconverged kernel attempt dcf09b87 is retained but excluded, with converged
replacement 320cdd53. No score-based filtering of the other completed methods.

## Completed result

Seven specialist configuration/seed runs completed, plus external efficacy
features, joint external graph training and full outer refits, width decisions,
and ten kernel variants. All model branches were audited. See
[specialist scope and seeds](2026-09-12-d6-specialist.md) and
[joint external training](2026-09-12-d6-joint-external.md).

The useful small signal is diversity of specialist representations. Their fixed
four-scope LGBM average scores 0.828992 alone; a prespecified 50% addition to the
previous global reference gives 0.822232 versus 0.825436. Repeating the D6-auxiliary
model across three seeds gives standalone 0.836105 versus individual 0.840989–0.844236.

Dedicated external efficacy features do not help. Joint external training improves
the direct graph output at public weight 0.30 (0.893099 to 0.884174), but its latent
LGBM does not beat the no-public refit (0.839558 versus 0.838061). An earlier positive
latent result using inner checkpoints disappears after full outer refitting.
Neither interval-width calibration nor count-fingerprint kernels materially helps.

The expanded 174-candidate Caruana pool is also negative: weight-crossfit D6
0.825436 to 0.825952, macro 0.582816 to 0.582945. Other endpoints are numerically
unchanged. Preserve the fixed-mean development candidate, but do not claim a
resolved improvement: its conditional group-bootstrap interval includes zero.
See [aggregate results](2026-09-12-d6-global-ensemble.md).

Historical receipts show improving D6 LB rank correlation without a comparable
absolute-error gain, but the three non-affine predictions cannot identify the
blind population moments used in the community's placement approach. See
[receipt-based diagnostic](2026-09-12-d6-submission-history.md). No blind moments
were invented, and no new submission was made.
