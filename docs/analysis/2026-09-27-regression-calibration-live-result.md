# Combined regression calibration Live result

**2026-09-27 (JST).** The combined 1A2, 2D6, and 3A4 affine calibration probe
was uploaded as prepared regression submission
`20260926T140601457029Z-bd3a1e40`, attempt
`20260927T051423190683Z-117605b9`. The official Live boards now show the
N283T row dated 2026-09-27 05:14 UTC. Audit run
`regression/20260927T053105312608Z-d49e8a6a` verifies the sent payload
SHA-256, all 750 molecule identities, the exact changed targets, and the
pre/post rows for all four endpoint boards and macro. The prior row is the
2026-09-26 10:00 UTC robust-D6 submission.

| Board | Prior ST-RAE | New ST-RAE | Change | Rank |
| --- | ---: | ---: | ---: | ---: |
| CYP1A2 | .508530 | .512345 | +.003815 | 40→42 |
| CYP2C9 | .359358 | .359358 | .000000 | 11→11 |
| CYP2D6 | .549341 | .549788 | +.000447 | 22→25 |
| CYP3A4 | .431253 | .431687 | +.000435 | 36→36 |
| Macro | .462120 | .463295 | +.001174 | 24→24 |

Lower ST-RAE is better. The macro increase equals the mean of the four
endpoint changes. The unchanged 2C9 payload and its scores match exactly.
All three calibration changes worsen Live ST-RAE, although D6 MAE improves
`.998980→.995539` and R² `.375395→.377681`. 1A2 and 3A4 MAE and R² worsen.
Rank correlations are unchanged on every board, as expected from the
rank-preserving positive affine maps. Ranks also depend on other entrants and
are secondary to the paired scores.

**Decision:** reject the three affine changes as improvements to the current
Live vector. Use the prior 2026-09-26 submission as the regression performance
control; preserve the new result as a calibration sensitivity observation.
The 1A2 and 3A4 OOF improvements did not transfer to the blind scored cohort.
The OOF screens and conditional bootstrap intervals were developmental and
did not account for calibration selection or the test-set activity shift.
The exact reason for the transfer gap is not identifiable from aggregate
leaderboard metrics alone. Do not resend either payload to recover the score.
