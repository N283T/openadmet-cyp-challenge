# Full-data Caruana deployment component

The retained original 24-candidate Caruana component was refitted for blinded test inference. This is one of the five equally weighted components of the selected OOF ensemble, not a newly evaluated OOF model.

- Source inner selections: regression `20260912T054406239560Z-96fce32b`.
- Full-data refit: regression `20260912T074545177210Z-0d561e17`.
- Independent membership and selection-transfer audit: regression `20260912T074646849858Z-1a37f7b6`.
- Auxiliary features: shared `20260912T074339967571Z-ff194fcf`; recipient-excluded cross-fitted train predictions and an all-measured-data auxiliary model for test predictions. Independent auxiliary audit: shared `20260912T074551656789Z-beecf03c`.

For each endpoint, the deployment weight of each candidate is the arithmetic mean of its five saved inner-selected fold weights, including zero when absent. Each positively weighted candidate is refitted once on **every available endpoint label**. No new candidate selection, weight fitting, calibration, or test-dependent tuning was performed. The original candidate parameters are retained, with a common deployment seed of 42 replacing the outer-fold seed offset. The fixed candidates have deterministic full-feature/full-row LightGBM fitting.

The eight original feature views remain Morgan2048, four predicted log2FC values, external classification embeddings600, external mixed-objective embeddings600, external tree scores20, and each of the latter three with the four log2FC values. Each view retains the original L1, quantile0.35, and quantile0.65 recipes: 300 trees, 15 leaves, minimum20 samples, learning rate0.05, L2=1. Only positive-weight candidates are materialized.

| Endpoint | Full training labels | Saved full-data models | Blinded test predictions |
|---|---:|---:|---:|
| CYP1A2 | 1412 | 20 | 750 |
| CYP2C9 | 1285 | 18 | 750 |
| CYP2D6 | 1493 | 19 | 750 |
| CYP3A4 | 2335 | 20 | 750 |

Independent validation confirms exact equality of all6,525 source label keys/values and training memberships, zero blinded-test overlap, exactly3,000 finite keyed predictions, and exact transfer of saved inner-selected weights. Model files and endpoint selection manifests are saved under each endpoint output directory. The keyed component prediction is `outputs/test-predictions.parquet`.

Reproduction: `pixi run python scripts/launch_full_caruana_refit.py --aux <full-aux.parquet>`, followed by `pixi run python scripts/launch_full_caruana_audit.py --run <new-run-id>`. The existing immutable configuration paths must be versioned before a new launch; completed runs remain intact.
