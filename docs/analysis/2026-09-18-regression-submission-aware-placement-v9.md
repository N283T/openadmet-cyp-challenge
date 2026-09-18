# Regression submission-aware placement v9

## Scored evidence

The September 12 unplaced submission and September 14 source-placement
submission establish that test-time distribution correction is material rather
than cosmetic. CYP2D6 moved from mean/SD `4.638/.405` to `3.570/.900`, and its
leaderboard ST-RAE improved from `1.334970` to `0.565520`. The September 17 flat
bank retained exactly those moments but changed ordering; D6 worsened to
`.582788` while Spearman fell from `.431249` to `.408851`. Thus the latest D6
loss is attributable primarily to ranking, not evidence that the validated
center should be discarded.

The other endpoint source placements remain useful. At the same fixed moments,
the September 17 model change moves 1A2 `.526455 -> .529845`, 2C9
`.366939 -> .359358`, and 3A4 `.444983 -> .431253`. This agrees with the
direction of their rank-metric changes and does not motivate a new location or
scale probe for those targets.

## D6 bank-delta audit

Run `0891022e` compares the old September 17 D6 simplex prediction with the v8
Caruana prediction using
`old + alpha * (new - old)`. Pooled OOF ST-RAE is best at the unamplified
`alpha=1`: `.814706`. Values `.5/1.25/1.5/2.0` score
`.815297/.814796/.815059/.816444`; `alpha=3` worsens to `.822369`.
Amplifying the small test-vector change merely to produce a larger submission
change is therefore rejected.

The v8 D6 OOF Pearson correlation is `.492241`. Multiplying it by the archived
leaderboard-inferred blind-label SD `1.599` gives a correlation-scaled prediction
SD of `.787094`. The adventurous v9 correction rounds this to the prespecified
`.80` sensitivity point while retaining the directly leaderboard-validated
ST-RAE center `3.57`. The conservative `.90` version is also preserved in the
run outputs and is numerically the v8 candidate.

## Candidate

Candidate `825d14db` uses source placements unchanged for 1A2, 2C9 and 3A4 and
uses mean/SD `3.57/.80` for D6. Relative to the September 17 submission, mean
absolute changes are `.09957/.06687/.08898/.07270`, and Spearman correlations
are `.99352/.99547/.99193/.99562` for 1A2/2C9/D6/3A4. Relative to the v8
candidate alone, the D6 SD update changes predictions by mean absolute `.07456`
and maximum `.36357`, without changing ranks. The earlier v8 D6 change was only
`.04704` versus the September 17 submission, so the new total D6 change is
materially larger without extrapolating the model delta.

This is an adventurous submission candidate, not a measured leaderboard
improvement. The mean `3.57` has direct scored support; the narrower `.80` spread
is a model-confidence adaptation derived from OOF correlation and public inferred
blind spread. No test label is used and no upload is performed.

## Selection status

At `2026-09-18 13:59 JST`, the user selected adaptive-D6 v9 as the next
regression submission candidate. Candidate `825d14db` remains unsubmitted during
the current cooldown. The selection pointer is recorded in
`configs/regression/next-submission-candidate-v9.json`.
