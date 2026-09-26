# CYP2D6 calibration evidence from scored submissions

## Question and scope

Can our previous scored submissions identify a better CYP2D6 correction for the
[robust-bound candidate](2026-09-26-regression-d6-robust-primary-bound-followup.md)?
This audit uses only the D6 column. No new model or external target is trained.
The current submitted September 24 D6 vector is identical to the September 14
source-placement D6 vector, with blind prediction mean `3.57` and sample SD
`.90`. The new robust blend has mean `3.57` but sample SD `.8835`.

## What the submissions actually compare

Recorded inventory `regression/20260926T094132305350Z-cebb315e` aligns all 750
identities from eight submitted D6 prediction files and calculates their moments
and pairwise affine residuals. The failed first run
`regression/20260926T094118047174Z-7587f06d` stopped while concatenating mixed
float32/float64 payloads and produced no usable scientific output; the corrected
run converts stored predictions to float64. Historical Live scores below come
from the existing [submission-history](2026-09-12-d6-submission-history.md),
[September 14 model-family](2026-09-14-regression-d6-model-family-screen.md),
and [September 18 verdict](../research/2026-09.md#rl0138--submission-aware-v9-regression-upload-and-its-leaderboard-verdict)
records. The September 24 D6 score is also present in the frozen
[Live snapshot](../../artifacts/tracked/leaderboard-snapshot-20260926/csv/live/regression/CYP2D6_pIC50_direct_inhibition.csv).

| Submitted D6 column | Blind mean | Sample SD | Live ST-RAE | Relationship |
| --- | ---: | ---: | ---: | --- |
| Sep 10 baseline | 4.593 | .312 | 1.4162 | Different model |
| Sep 11 TabICL replacement | 4.629 | .317 | 1.4125 | Different model |
| Sep 12 ensemble | 4.638 | .405 | 1.3350 | Different model |
| Sep 14 placed control | 3.570 | .900 | **.5655** | Scored placement and model change |
| Sep 14 MoLFormer RF addition | 3.570 | .886 | .5725 | Model/ranks changed |
| Sep 17 flat bank | 3.570 | .900 | .5828 | Model/ranks changed |
| Sep 18 adaptive v9 | 3.570 | .800 | .5945 | Model/ranks also changed |
| Sep 24 endpoint-wise combination | 3.570 | .900 | **.5655** | D6 column exactly repeats Sep 14 |

Only the Sep 14 placed control and Sep 24 endpoint-wise combination are the
**same D6 vector** (maximum affine reconstruction error `2.7e-15`). All other
pairs change molecular ranks or have material non-affine residuals. For example,
Sep 17 to Sep 18 has Pearson `.9979`, but affine residual RMSE `.0521` and
maximum error `.179`; Sep 12 to Sep 14 has RMSE `.1503`. Thus none of the
different-scoring pairs isolates a correction of one fixed model. Treating the
Sep 17→18 score loss as an exact `.90→.80` calibration trial would conflate
spread with the simultaneous prediction changes. Multiple submissions do not
identify an optimal new D6 mean or SD from our own receipts.

The historical data do support preserving the known competitive `3.57/.90`
placement as a **conservative reference**. That is evidence about a previously
scored vector, not proof that the same moments are optimal for the new model.
The Live and one-time full-test Interim boards also score different molecule
populations; do not pool their scores as repeated calibration measurements.

## Label-free final-placement test for the new D6 blend

Run `regression/20260926T094321680608Z-a298391f` tests three fully specified
post-blend transformations on the new robust D6 OOF, applied separately within
each held fold without using its labels:

1. Keep the original blend.
2. Restore its mean and SD to the matched control's prediction moments.
3. Match the entire empirical control-prediction distribution by quantile.

| Final placement | Parent ST-RAE | Scaffold ST-RAE |
| --- | ---: | ---: |
| Original blend | **.819829** | .810041 |
| Restore control SD | .820809 | .809318 |
| Match control quantiles | .823512 | **.807774** |

Both distribution-restoring changes worsen parent and improve scaffold.
Their parent/scaffold disagreement and small changes do not support a new
primary calibration rule. These OOF transformations also cannot validate the
blind D6 mean `3.57`; the labelled D6 OOF population has a different activity
mix from the challenge test.

Run `regression/20260926T094437902447Z-a8262487` nevertheless packages the
minimal blind sensitivity variant that restores the **final** robust blend's
sample SD from `.88348` to `.90000` while keeping mean `3.57`. It changes the
750 D6 predictions by mean absolute `.0124`, maximum `.0644` pIC50, preserves
their ranks, and leaves other endpoints unchanged. It is a reviewable
calibration alternative, **not a demonstrated improvement**. No live upload
was made.

## Decision

Keep the original robust-bound candidate as the primary exploratory model.
The historical receipts do not justify fitting a more elaborate calibration
curve or changing the `3.57` mean. The exact-SD-`.90` sensitivity variant is
available if a future controlled submission tests final spread. A genuinely
informative calibration comparison would hold the 750 D6 prediction ranks
fixed and submit prespecified location/scale variants, subject to leaderboard
limits; model-changing historical submissions cannot answer that question.
Compact outputs and the sensitivity payload are in
[`artifacts/tracked/regression-d6-calibration-history-20260926/`](../../artifacts/tracked/regression-d6-calibration-history-20260926/README.md).
