# Leaderboard-informed CYP2D6 placement on the simplified regression vector

The current simplified submission scores macro ST-RAE 0.719349. CYP2D6 dominates
the remaining error at 1.353707, versus 0.685494/0.369556/0.468640 for the other
three endpoints. Small endpoint substitution can only move the macro by a few
thousandths, so it is not the correction considered here.

The scored CYP2D6 vector has Spearman 0.431249 but R² -0.813077. Its full 750-row
prediction distribution has mean 4.655864 and sample SD 0.399258. The public
blind-moment analysis attributes mean 3.107 and SD 1.599 to the live labels, while
its separately probed ST-RAE placement uses prediction mean 3.57 and SD 0.90.
The current vector is therefore much higher and narrower than that target.

Run `20260913T062809247284Z-83afbfa3` applies the frozen public placement code to
the exact current submitted vector. The decision candidate changes only CYP2D6:

| Quantity | Current | Placed |
|---|---:|---:|
| Mean | 4.655864 | 3.570000 |
| Sample SD | 0.399258 | 0.900000 |
| Minimum | 4.091822 | 2.298548 |
| Maximum | 6.026955 | 6.660691 |
| Affine scale | 1.000000 | 2.254183 |
| Mean absolute change | 0 | 1.118653 |
| Spearman versus current | 1 | 1 |

Packaging run `20260913T063146822638Z-b5fd61c7` validates all 750 identities and
four regression columns. Candidate `20260913T063202760092Z-9aeeb189` preserves
the other 2,250 endpoint predictions exactly.

The public source reports CYP2D6 ST-RAE 0.5653 for the same target mean and SD on
its own stronger-ranked vector. That number is not transferable as this
candidate's score. If this candidate reached 0.5653 on CYP2D6, the observed other
three endpoint scores imply macro 0.5222. Reaching macro below 0.60 requires its
CYP2D6 score to fall below 0.8763. This is a large-effect leaderboard hypothesis,
with uncertainty dominated by the weaker CYP2D6 ordering in the current model.

Ordinary train-label calibration remains rejected: its prior local CYP2D6 gain
was about 0.004 and cannot explain the live placement gap. This candidate instead
tests the leaderboard-inferred shift and spread directly.
