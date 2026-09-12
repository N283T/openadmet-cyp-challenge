# Jeremy FT embedding + log2fc Rustuna readout

Training run `20260912T072056445264Z-e6b6cbcb` evaluates 24 LightGBM trials per outer fold and endpoint. The grouped 20% inner holdout matches the existing Caruana membership exactly; anchor groups remain in inner training. The saved winner is the minimum official inner ST-RAE in all 20 fits, and all 480 trial scores are finite. Winners refit on all outer training labels.

The FT encoder previously used labels from the readout's inner holdout, so inner HPO scores can be optimistic. This is not nested encoder selection. Each encoder's outer holdout was verified against parent assignments; those outer labels were excluded from encoder fitting, checkpoint selection and readout HPO. Upstream pretraining still lacks a molecule manifest, so this remains an upstream-model diagnostic.

The comparison retains the fixed LGBM readout and uses the same prescribed equal three-way mean: original 24-candidate Caruana + frozen Jeremy TabPFN + FT/log2fc LGBM. No ensemble weights were fitted to OOF predictions.

Official standalone audit: `20260912T072848133969Z-be4d2d35`. Inner-selection verification and fixed-blend comparison: `20260912T072850188008Z-da7355fc`.

Tuning improves overall and neighbor scores for both the standalone readout and the prescribed three-model blend. Improvements come from 2C9 and 3A4; 1A2 and 2D6 worsen slightly. The three-model result is 0.603423, so it does not itself achieve the sub-0.60 target. The separately evaluated fixed five-model ensemble is outside this comparison.

| Method | All ST-RAE | Neighbors ST-RAE |
|---|---:|---:|
| fixed_lgbm | 0.628743 | 0.747805 |
| tuned_lgbm | 0.624940 | 0.713751 |
| fixed_three | 0.605039 | 0.696275 |
| tuned_three | 0.603423 | 0.687040 |

| Endpoint | Fixed LGBM | Tuned LGBM | Fixed three | Tuned three |
|---|---:|---:|---:|---:|
| CYP1A2 | 0.711147 | 0.720947 | 0.686811 | 0.689106 |
| CYP2C9 | 0.505252 | 0.494596 | 0.479764 | 0.475299 |
| CYP2D6 | 0.851538 | 0.853844 | 0.838174 | 0.839918 |
| CYP3A4 | 0.447035 | 0.430372 | 0.415407 | 0.409371 |
