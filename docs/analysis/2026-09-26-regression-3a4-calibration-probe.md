# CYP3A4 calibration-only probe

**2026-09-26 (JST).** Following the classification-transfer model search,
screen a small affine correction of the *scored* CYP3A4 regression vector.
The model and molecule ranking stay fixed. The current scored 3A4 column has
Live ST-RAE `.431253`, Interim ST-RAE `.3868`, and blind mean/sample SD
`4.880000/1.082472`. Its location and spread originally came from the source
author's leaderboard-informed placement constants; they are not established
optimal moments for our scored model or its 750-row test vector.

## OOF screen

Run `regression/20260926T140255759268Z-7ddd1205` applies fixed shifts
`−.2…+.2` and scales `.85…1.15` around the prediction mean of the exact
September 17 scored 3A4 bank on 2,300 parent-fold OOF rows. Its identity
ST-RAE is `.394567`. The grid minimum is shift `+.10`, scale `.925`:
ST-RAE `.386083` (change `−.008484`), but MAE worsens
`.467143→.492401`. A shift of `+.05` without narrowing scores `.389082`.
Spearman is unchanged by these monotone maps.

Run `regression/20260926T140426997076Z-e44344d9` evaluates the grid
winner and two smaller fixed corrections by fold and by 1,000 bootstrap
replicates sampling 1,083 molecular-series groups:

| Shift / scale | OOF ST-RAE change | Fold wins | Conditional 95% series interval |
| --- | ---: | ---: | ---: |
| `+.05 / 1.0` | −.005485 | 5/5 | [−.009153, −.002031] |
| `+.10 / .925` | −.008484 | 4/5 | [−.014793, −.002644] |
| `+.05 / .9625` | −.006659 | 5/5 | [−.009835, −.003809] |

These bootstrap intervals condition on the saved OOF predictions and the
selected variants. They do not correct for the grid search, refit the bank,
or measure transfer to blind data. OOF prediction mean/sample SD are
`4.154160/.805189`, unlike the submitted `4.88/1.082472`; test selection
targets 3A4 active analogs and can change its optimal calibration. The
current source placement was itself informed by prior leaderboard probes.

## Prepared candidate

Choose the modest `+.05 / .9625` map for a prospective sensitivity test: it
improves all five folds and changes the established blind placement less than
the screen winner. Run `regression/20260926T140539644749Z-faeac48c` applies
it around the scored blind mean. The 750-row 3A4 mean/sample SD become
`4.930000/1.041879`; mean absolute prediction change is `.052227` pIC50,
maximum `.168446`. The ranking, molecular identities, and model are fixed.
The run validates a 3A4-only file and a combined file preserving the already
prepared 1A2 narrowing and D6 SD-`.90` probes. CYP2C9 stays identical.

3A4-only candidate `20260926T140553016103Z-e1e32ffc` and combined candidate
`20260926T140557971440Z-00366073` are recorded. The combined candidate was
prepared as submission `20260926T140601457029Z-bd3a1e40`, SHA-256
`ea5151a2be100576a70cda36ae600ebb8871aff565cab93789236879307dc20d`;
the dry run passed 750-identity, schema, and profile checks. **No upload was
made.** Target-specific leaderboard results can evaluate the three corrections
separately when the submission cooldown ends, using the September 26 scored
submission as the common control. OOF improvement is evidence for trying the
probe, not a claimed Live improvement.
