# Adoption of the Rustuna primary-screen model

The user explicitly requested replacing the auxiliary model with the tuned recipe.
This supersedes the previous experiment's decision to preserve it only as a candidate.
The current regression feature source is Rustuna trial 14 from
`20260911T120540987924Z-3f5823d1`. Its settings are fixed during completion of all
ten nested fold-pair models. The existing four fold-zero pairs are reused; the
six remaining pairs are trained once with the same settings.

The canonical registry is `configs/regression/current-log2fc.json`, consumed by
`scripts/launch_parent_readout.py`. It carries the explicit selected recipe, split,
ten model run IDs and a validated package run. Package validation verifies selected
hyperparameters, pair coverage, held-out membership, complete molecule IDs and
finite four-head predictions. It produces all ten nested prediction tables,
four-model OOF predictions for 6,146 non-test compounds and ten-model averaged
auxiliary predictions for the 750 test compounds.

The four remaining LightGBM readout folds are also updated. The existing tuned
fold-zero readout is reused. `configs/regression/current-parent-readout.json` points
to these runs and is consumed by `scripts/launch_tabicl_parent.py`. This updates
the source for future TabICL experiments; it does not imply that new TabICL fits
have been run. Historical configs, model assets and prior submissions remain intact.

Fold 0 was excluded from hyperparameter selection. The other folds were involved
in that selection, so the completed five-fold metrics describe the adopted recipe
after selection and are not a fresh nested-HPO estimate. Adoption follows the
user's decision rather than another search over these results.

## Completed assets and validation

The six additional auxiliary runs (pairs 1/2, 1/3, 1/4, 2/3, 2/4, 3/4) are:

- `20260911T124637213478Z-68a164d5`
- `20260911T124755027333Z-d1af8608`
- `20260911T124921417170Z-9770d064`
- `20260911T125041415499Z-d0d0f393`
- `20260911T125207635438Z-89264e78`
- `20260911T125345536596Z-0cd216fa`

Package `20260911T125618165690Z-4e52f1df` validates all ten pairs, 6,896 finite
prediction rows per model, 6,146 OOF compounds and 750 test compounds. Initial
packaging attempt `20260911T125539789954Z-0cbd4373` failed while writing the registry
because a local dataframe shadowed the selected-configuration variable. It is
retained; the fix changes packaging only and requires no model retraining.

Updated readout runs for folds 1–4:

- `20260911T125630683489Z-e2f0a4b2`
- `20260911T125723862156Z-08fd7a35`
- `20260911T125818040377Z-e0d35c62`
- `20260911T125911665760Z-06163bb5`

Fold 0 reuses `20260911T124052504276Z-c22b0df7`. Audit
`20260911T130021181456Z-60011f3a` validates the full paired OOF set and registers
576 metrics. Full-feature parent-retained results:

| Subset | Previous ST-RAE | Adopted ST-RAE | Previous MAE | Adopted MAE |
|---|---:|---:|---:|---:|
| All nonanchors | 0.687660 | 0.684192 | 0.537846 | 0.537255 |
| Hit neighbors | 0.836776 | 0.850116 | 0.547848 | 0.557958 |

These metrics do not change the user's adoption decision. The reference registry,
README and launchers now consistently select the adopted auxiliary models and
matching readout inputs. Existing UMAP experiment launchers remain historical.
Checks and 44 DB-enabled tests passed. No submission was made.
