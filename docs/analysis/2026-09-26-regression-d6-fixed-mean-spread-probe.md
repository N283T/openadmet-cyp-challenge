# CYP2D6 next submission: fixed-mean spread probe

## Rationale and payload

The submitted robust-bound D6 blend improved the [Live D6 ST-RAE](2026-09-26-regression-d6-robust-live-result.md)
from `.565520` to `.549341`. Its blind prediction mean is `3.57` and sample SD
`.8834767`. Historical submissions do not isolate a calibration optimum because
their D6 model predictions also changed. A controlled next submission can hold
all 750 D6 ranks and the mean fixed, changing only the sample SD to `.90`, the
previously scored control spread. The multiplier around mean `3.57` is
`1.0187025685`.

Run `regression/20260926T094437902447Z-a8262487` created the variant before
the robust blend was uploaded. Its source candidate SHA-256
`399cd07c66c0d4e7e1fb468b6a36a471d1572951e7f10e324e0b597d216c7b29`
exactly matches the subsequently scored submission. Candidate
`20260926T120120315669Z-23ae8720` and prepared regression submission
`20260926T120124524160Z-01cea54f` package the variant; the submission
dry run passes 750-row, schema, hash and profile validation. Payload SHA-256 is
`25074ee65b258b24eda2dabf268b12eabb91a813a5abeb65f605735dcbc400fc`.
The mean absolute D6 change is `.012424` pIC50 and maximum change `.064360`.
The other three regression columns are byte-identical at the value level to
the scored payload. This prepared submission has **not** been uploaded.

## Matched OOF sensitivity

The earlier [calibration audit](2026-09-26-regression-d6-submission-calibration-audit.md)
restored both the control fold mean and SD. That does not exactly match the
proposed blind transformation, which preserves the robust blend's mean. Run
`regression/20260926T120245745572Z-ea66533b` applies the exact blind multiplier
around each held fold's own prediction mean, without using held labels. It
checks the two blind payload hashes, exact test identities, unchanged other
endpoints and the proposed affine D6 relationship. Full keyed OOF predictions
are retained in the run; [compact metrics](../../artifacts/tracked/regression-d6-fixed-mean-spread-20260926/README.md)
are tracked.

| Split | Submitted ST-RAE | SD .90 ST-RAE | Change | Fold direction |
| --- | ---: | ---: | ---: | --- |
| Parent | .819829 | .820562 | +.000733 | Worse in 5/5 |
| Scaffold | .810041 | .809522 | −.000520 | Better in 4/5 |

The paired OOF result is small and split-dependent. On the blind 750 rows the
transformation is strictly increasing, so their D6 order is unchanged. The
OOF Spearman differs slightly because independent fold-wise centering can
change order between folds. The labelled OOF activity mix also differs from
the blind test, and neither split estimates the Live ST-RAE response reliably
enough to call `.90` an upgrade.

## Decision

Keep the currently scored robust blend as the performance control. The `.90`
variant is a **calibration-only prospective probe**: if uploaded, its paired
Live score would reveal the local direction of spread change for this fixed
D6 vector. It should not be described as an offline-proven improvement. No
mean shift, nonlinear map or new model is bundled into this probe, so the
result remains interpretable.
