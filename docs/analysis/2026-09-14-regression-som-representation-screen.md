# Site-of-metabolism representations for direct-inhibition regression

## Question and corrected design

The classification SoM screen produced reusable SMARTCyp, challenge-disjoint
ATTNSOM, and frozen CYP-MAP features for all challenge compounds. This follow-up
asks whether those representations improve direct-inhibition pIC50 regression.

The first registered attempt, `20260914T111246035102Z-cfe292cc`, incorrectly used
the legacy shared UMAP/scaffold comparison and Mordred+Morgan control. It was
interrupted before producing a result and is excluded. The valid paired run is
`20260914T111634758781Z-477d7e4e`. It uses the retained parent split
`20260911T104703655337Z-a3923255`, in which hit-seed anchors always remain in
training, and exactly reproduces the existing frozen-600 plus nested-log2fc
LightGBM control from `20260912T083010595992Z-1e7001de`.

Every arm uses the same four endpoints, outer folds, labels, confidence bounds,
base features, LightGBM parameters, and validation membership. The control is
bitwise identical over all 6,384 keyed OOF rows. The additions are frozen and do
not use direct-inhibition labels during representation training. CYP-MAP has 155
connectivity overlaps with challenge training compounds but none with the blind
test; its source labels are metabolism tasks rather than direct-inhibition pIC50.

## Paired readout result

Lower ST-RAE is better.

| Added block | 1A2 | 2C9 | 2D6 | 3A4 | Macro |
|---|---:|---:|---:|---:|---:|
| None (control) | .691998 | .494360 | .842287 | .442934 | .617895 |
| SMARTCyp | .688915 | **.489679** | **.842062** | .436766 | **.614356** |
| ATTNSOM | .691753 | .499572 | .846795 | .444762 | .620721 |
| CYP-MAP embedding | .696497 | .493571 | .847400 | .437795 | .618816 |
| CYP-MAP all | .696407 | .493488 | .848612 | **.435728** | .618559 |
| SMARTCyp + ATTNSOM | **.688465** | .492015 | .846646 | .442932 | .617514 |
| All embedding-oriented blocks | .689459 | .496129 | .846696 | .437634 | .617479 |
| All blocks | .701734 | .495741 | .849469 | .441241 | .622046 |

SMARTCyp is the only clean macro improvement and lowers ST-RAE on all four
endpoints. Its 2D6 change is negligible (`-.000225`); most of the improvement
comes from 3A4 and 2C9. ATTNSOM and CYP-MAP do not provide a general direct-
inhibition gain when appended at full feature weight. The high-dimensional union
is actively harmful.

## Addition to retained ensembles

Run `20260914T112745470705Z-1c0d6aa1` verifies all keys, labels, confidence bounds,
and parent memberships, then applies fixed 10%, 25%, and 50% blends to the latest
development OOF vector (`20260912T083627802028Z-ae562a59`) and the submitted
five-component reference.

No mechanism representation improves the latest macro score. The smallest
SMARTCyp addition changes `.586931` to `.587254` at 10%; larger weights are worse.
Endpoint-specific point changes are too small and inconsistent for deployment:
10% SMARTCyp changes 2D6 `.822681` to `.822316`, while 25% all-embedding features
change 1A2 `.670809` to `.669357`. The same additions worsen 2C9 and 3A4. CYP-MAP
embedding improves the small neighbor-only macro at some weights, but this does
not survive the all-row objective and is not a candidate-selection result.

## Decision

Retain SMARTCyp as a useful low-dimensional regression feature and evidence that
SoM-related accessibility correlates with direct inhibition. Do not add any of
these readouts to the current regression candidate and do not create a test
inference package or submission. The experiment does not supply the genuinely
new CYP2D6 ordering signal required after the successful leaderboard placement.

