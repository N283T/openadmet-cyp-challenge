# Unsubmitted D6 AP-trunk probe: blind call review

The user authorized inspection only and explicitly prohibited submission.

Selected full-development recipe: `['core_separate_fe--trunk--lgbm_binary']`; trunk rank weight 0.50.
Core-pocket, Fe and remaining HEM pair embeddings are aggregated separately and read with LightGBM. Five raw fold predictions are averaged, then ranked across the 750 blind compounds. The final score blends this rank with the existing D6 ensemble rank.

D6: 49 positives; 28 retained, 21 added, 21 removed.
A4: 275 positives, unchanged from the last sent 275-call probe. This does not restore the earlier best-performing 249-call vector.

The previously reported AP .296860→.318808 is the conditional held-fold procedure result, not a validation score measured for this exact final blind vector. The full-development inner AP is a recipe-selection diagnostic, not independent confirmation.

## Changed blind calls

| Change | Molecule | Control position | Candidate position | Weighted positive-neighbor fraction | Trunk-fold positive votes |
|---|---|---:|---:|---:|---:|
| added | OCNT-2534965 | 58 | 16 | 0.432 | 5/5 |
| added | OCNT-2535110 | 67 | 17 | 0.000 | 5/5 |
| added | OCNT-2534916 | 64 | 18 | 0.415 | 4/5 |
| added | OCNT-2535286 | 92 | 20 | 0.333 | 5/5 |
| added | OCNT-2535563 | 96 | 21 | 0.105 | 5/5 |
| added | OCNT-2534939 | 88 | 22 | 0.461 | 4/5 |
| added | OCNT-2534912 | 59 | 23 | 0.533 | 4/5 |
| added | OCNT-2535053 | 73 | 25 | 0.091 | 2/5 |
| added | OCNT-2535810 | 101 | 30 | 0.232 | 5/5 |
| added | OCNT-2534870 | 66 | 31 | 0.468 | 4/5 |
| added | OCNT-2535348 | 100 | 32 | 0.487 | 4/5 |
| added | OCNT-2535538 | 84 | 33 | 0.350 | 4/5 |
| added | OCNT-2535287 | 68 | 34 | 0.420 | 3/5 |
| added | OCNT-2534938 | 115 | 37 | 0.416 | 4/5 |
| added | OCNT-2534956 | 121 | 38 | 0.204 | 3/5 |
| added | OCNT-2535643 | 95 | 39 | 0.434 | 3/5 |
| added | OCNT-2535598 | 125 | 42 | 0.199 | 2/5 |
| added | OCNT-2534921 | 134 | 43 | 0.339 | 4/5 |
| added | OCNT-2535422 | 99 | 44 | 0.096 | 3/5 |
| added | OCNT-2534943 | 130 | 45 | 0.430 | 3/5 |
| added | OCNT-2535667 | 86 | 49 | 0.239 | 2/5 |
| removed | OCNT-2535279 | 15 | 52 | 0.486 | 2/5 |
| removed | OCNT-2535536 | 39 | 60 | 0.177 | 2/5 |
| removed | OCNT-2535435 | 7 | 64 | 0.485 | 2/5 |
| removed | OCNT-2535151 | 49 | 66 | 0.453 | 2/5 |
| removed | OCNT-2535147 | 32 | 68 | 0.177 | 3/5 |
| removed | OCNT-2535543 | 28 | 70 | 0.093 | 1/5 |
| removed | OCNT-2535496 | 1 | 92 | 0.677 | 1/5 |
| removed | OCNT-0493952 | 36 | 97 | 0.293 | 1/5 |
| removed | OCNT-2535495 | 35 | 102 | 0.445 | 1/5 |
| removed | OCNT-2535280 | 5 | 121 | 0.595 | 0/5 |
| removed | OCNT-2535178 | 27 | 143 | 0.428 | 0/5 |
| removed | OCNT-2535498 | 9 | 152 | 0.467 | 0/5 |
| removed | OCNT-2535126 | 41 | 171 | 0.510 | 0/5 |
| removed | OCNT-2535497 | 22 | 184 | 0.383 | 0/5 |
| removed | OCNT-2535278 | 8 | 190 | 0.486 | 0/5 |
| removed | OCNT-2535240 | 23 | 224 | 0.709 | 0/5 |
| removed | OCNT-2535520 | 40 | 240 | 0.605 | 0/5 |
| removed | OCNT-2535611 | 44 | 253 | 0.414 | 0/5 |
| removed | OCNT-2534924 | 37 | 264 | 0.079 | 0/5 |
| removed | OCNT-2535368 | 47 | 265 | 0.326 | 0/5 |
| removed | OCNT-2534824 | 24 | 273 | 0.269 | 0/5 |

Full keyed scores, SMILES, nearest positive/negative similarities, weaker-weight sensitivity and trunk-fold sensitivity are retained alongside this report.

No upload or model training was performed. No full-pool features were added.
