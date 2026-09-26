# CYP2D6 robust-bound submission: Live leaderboard result

## Submission and audit

The user authorized submitting the prepared 50:50 CYP2D6 robust-bound blend
with blind prediction mean `3.57` and sample SD `.88348`. Candidate
`20260926T090542937263Z-099c328f` was prepared as regression submission
`20260926T090548323924Z-46188564`, payload SHA-256
`399cd07c66c0d4e7e1fb468b6a36a471d1572951e7f10e324e0b597d216c7b29`.
Attempt `20260926T100008522917Z-0ab7c954` received the official service
receipt for N283T on the Regression Prediction track. The receipt says the
submission was processing, not that it had scored.

The reflected N283T row dated **2026-09-26 10:00 UTC** appeared in the Live D6
leaderboard snapshot `20260926T100306712074Z-6e0eb8fc`; the Live macro row
appeared in `20260926T100322667181Z-7a237ebe`. Submission event
`ba93220c876c47e29cfb3f5611b78a7e` records the accepted D6 snapshot.
The comparison uses pre-upload Live snapshots taken automatically in the same
attempt, not a prior board from a different scoring population.

Recorded audit `regression/20260926T100510447246Z-e2a865c9` checks the receipt,
payload hash, 750 molecule identities and exact preservation of the other
three regression columns. It matches the pre/post N283T rows by their actual
submission timestamps and confirms the macro ST-RAE change is exactly one
quarter of the D6 ST-RAE change. Compact
[audit metrics](../../artifacts/tracked/regression-d6-robust-live-20260926/README.md)
are tracked. Initial audit run `regression/20260926T100435921054Z-60a00a3d`
failed because the macro leaderboard prefixes metric columns with `MA-`;
it yielded no usable scientific result. The corrected recorded run above
completed successfully.

| Live metric | September 24 control | September 26 robust blend | Change |
| --- | ---: | ---: | ---: |
| D6 ST-RAE ↓ | .565520 | **.549341** | **−.016179** |
| D6 MAE ↓ | 1.009140 | **.998980** | −.010160 |
| D6 R² ↑ | .364216 | **.375395** | +.011179 |
| D6 Spearman ↑ | .431249 | **.441151** | +.009902 |
| Macro ST-RAE ↓ | .466165 | **.462120** | −.004045 |

The D6 Live rank moved `31→22` and macro rank `27→25`; ranks depend on other
entrants and are secondary to the paired metric change. This is a real
prospective improvement on the **Live scored half** of the blind test. The
Interim full-test board was frozen before this submission, so it does not
evaluate this candidate. Aggregate leaderboard metrics cannot identify which
compounds improved or whether the same gain will hold on the hidden half.

## Calibration decision

The submitted final SD was `.88348`, slightly below the historically scored
`.90` reference. The separate
[calibration audit](2026-09-26-regression-d6-submission-calibration-audit.md)
created an exact-SD-`.90` variant, but its parent/scaffold OOF results disagree
on direction. The Live result supports the **submitted combination of model and
calibration**, not a claim that `.88348` is the optimal spread. Do not alter
the newly scored D6 control solely to force `.90`; any spread comparison must
hold the D6 prediction ranks fixed in another prospective submission. No such
second submission was made in this turn.
