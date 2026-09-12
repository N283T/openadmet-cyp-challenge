# Source-author placement: assumptions and test scenarios

The author's placement is a substantial change to predicted level and spread.
Applied to our already submitted five-component vector, its D6 ST-RAE setting
changes mean/SD from **4.6382/0.4054 to3.57/0.90**, without changing ranks. A
separate D6-only scenario preserves every other endpoint value exactly.
These are prepared prediction scenarios, not measured improvements; nothing was
uploaded by these runs.

## Formula and public evidence

Pinned Workbench commit `a1d1b46e290e62014ee3f0dfde5667eb5e93282b`, frozen in
`5ed4d83b`, contains
[cyp_recalibrate.py](https://github.com/SuperCowPowers/workbench/blob/a1d1b46e290e62014ee3f0dfde5667eb5e93282b/ml_pipelines/OpenADMET/cyp/scripts/cyp_recalibrate.py).
Its formula is:

`placed = max(1, target_mean + target_sd * (prediction - mean(prediction)) / sample_sd(prediction))`.

The centering and sample SD use all750 supplied test predictions. Its default
spread is `rho * reported_blind_sd`; `--strae` instead directly sets D6 SD to0.90
and mean to3.57. There is no maximum clipping. A floor can introduce ties.

| Endpoint | Author-reported blind mean | Author-reported blind SD | Author aux100 Pearson | Default target prediction SD |
|---|---:|---:|---:|---:|
| 1A2 | 4.412 | 1.553 | 0.753 | 1.169409 |
| 2C9 | 4.830 | 1.101 | 0.838 | 0.922638 |
| 2D6 | 3.107 | 1.599 | 0.678 | 1.084122 |
| 3A4 | 4.880 | 1.272 | 0.851 | 1.082472 |

The constants come from the author's scored affine submission probes. The D6
3.57/0.90 setting is explicitly sampled against leaderboard ST-RAE, not learned
from challenge training labels. `SOLVED_PEARSON` belongs to the author's aux100
model. An alternate source path multiplies each new model's OOF Pearson by the
author's ratios1.32/1.23/1.66/1.07 for1A2/2C9/2D6/3A4. Those ratios are
model/split observations, not established transport constants.

The script calls its blind-moment solve exact, but also says R² scores only a live
subset of750 compounds. Its implementation reads moments from full submission
vectors without selecting that scored subset. Exact R² moment equations require
moments and labels on the same cohort. Together with rounded board R² values,
this prevents independently treating the published constants as verified label
moments from source text alone. We retain them as clearly declared assumptions.
The source comment also says ST-RAE “fell” while its quoted sequence rises from
0.565 to0.694; the quoted values and surrounding prose describe deterioration.

## Executed scenarios

`scripts/scp_placement_scenarios.py` extracts the audited constants and original
`place` function from the frozen source AST. A minimal column adapter reproduces
Pandas mean/sample-SD behavior in float64. It applies the exact resulting affine
parameters and floor, saves all identities, and verifies monotonicity and exact
unchanged controls.

- Native ensemble: final run `20260912T113338205386Z-46754704`, based on
  `0c4a365a`. Includes source solved-R² and ST-RAE placements, our measured OOF
  Pearson times the author's ratios, D6-only correction,25/50/75% partial
  placements, and a declared D6 center3.32/3.57/3.82 × SD0.675/0.90/1.125 grid.
- Previously submitted vector: run `20260912T113337187016Z-d6e4f5d9`, based on
  package `45d23980`, with all3000 values independently matched exactly against
  submitted artifact `eaafa860` (sent attempt `3d2e5c50`). It includes raw,
  source-solved placements, D6-only and sensitivity scenarios. No mismatched
  native-model OOF Pearson is assigned to this older submitted vector.
- Initial native run `8e86e04f` remains recorded. The final run explicitly promotes
  predictions to float64 and preserves raw controls exactly rather than allowing
  float32 subtract/add round-off. It also adds the requested D6-only scenario.

| Vector / scenario | D6 mean | D6 sample SD | Scale | Offset | Floored |
|---|---:|---:|---:|---:|---:|
| Submitted raw | 4.638156 | 0.405434 | 1 | 0 | 0 |
| Submitted source ST-RAE / D6-only | 3.570000 | 0.900000 | 2.219846 | −1.068156 | 0 |
| Native raw | 4.687115 | 0.360996 | 1 | 0 | 0 |
| Native source ST-RAE / D6-only | 3.570000 | 0.900000 | 2.493104 | −1.117115 | 0 |

The corrected submitted D6 range is2.16949–6.55948, with mean absolute change
1.09916pIC50. Its ranks remain identical. All-endpoint source ST-RAE placement
floors one submitted1A2 prediction; the D6-only scenario preserves all2250 other
endpoint values exactly, including their tails.

Each run saves `scenario-test-predictions.parquet`, keyed by scenario, compound
and endpoint; `scenario-distributions.csv`; and `assumptions.json`. For a future
same-vector comparison, select `raw` and `d6_only_author_strae` from `d6e4f5d9`.
The all-endpoint alternative is `author_solved_strae`. These are not submission
records and carry no assumed or measured test-score column.
