# CYP2D6 submission-calibration audit

Compact evidence from the [audit](../../../docs/analysis/2026-09-26-regression-d6-submission-calibration-audit.md).
Full keyed predictions and immutable source snapshots remain in recorded runs
`regression/20260926T094132305350Z-cebb315e`,
`regression/20260926T094321680608Z-a298391f` and
`regression/20260926T094437902447Z-a8262487`.

| File | SHA-256 |
| --- | --- |
| `submission-moments.csv` | `f4027f109df670b4e5d7843e2cc3059485cf9b1d10d2fc2430a75cbc7898f37f` |
| `affine-comparisons.csv` | `d6314397cf97c990cbf1d4eb6cf3ff8a4b761d0a27c3b8c8887a8684ae6fc778` |
| `placement-metrics.csv` | `c4bd56a94720893fddabbd0255b3db84fbb79de85b3056f84ce0f645e4da691b` |
| `sd90-summary.json` | `c0d4df1b647d8d613d4f385e21902ff896c64439d10c32ff5e5284e4d6bd4127` |
| `candidate-sd90.csv` | `25074ee65b258b24eda2dabf268b12eabb91a813a5abeb65f605735dcbc400fc` |

`candidate-sd90.csv` is a sensitivity variant, not a selected submission.
It restores the final D6 blend to mean `3.57`, sample SD `.90`, and changes no
other endpoint. It has not been uploaded.
