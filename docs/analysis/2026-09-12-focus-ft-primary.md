# Replacing the primary features in the tuned FT component

Run `20260912T083421950497Z-036b9595` preserves the original Jeremy medium
challenge-fine-tuned embeddings from `20260912T060538153709Z-8d434eee` and all
20 previously selected per-fold/per-endpoint LGBM parameter sets from
`20260912T072056445264Z-e6b6cbcb`. Only the four nested primary-screen prediction
features change to the pair-excluded LGBM source
`20260912T082554377951Z-d272774c`. No parameters were retuned.

All 20 readouts are refit on their outer-training rows, with three CPU threads.
The comparison includes the submitted five-model ensemble and fixed 25%/50%
additions of the new component, plus a direct replacement of its old one-fifth
component. These weights were specified before evaluating this experiment.

| Model | CYP1A2 | CYP2C9 | CYP2D6 | CYP3A4 | Macro |
|---|---:|---:|---:|---:|---:|
| Old component | .720947 | .494596 | .853844 | .430372 | .624940 |
| New primary features, same parameters | .688647 | .467318 | .836937 | .414475 | .601844 |
| Submitted ensemble | .683335 | .463781 | .834741 | .403800 | .596414 |
| Replace old component one-for-one | .677241 | .458803 | .832048 | .399764 | .591964 |
| Add new component at 25% | .676253 | .453807 | .829758 | .396613 | .589108 |
| Add new component at 50% | .674383 | .450451 | .827964 | .395930 | .587182 |

Every endpoint improves when only the primary features are replaced. The
neighbor-subset macro also improves from .687208 for the submitted ensemble to
.676870 at 25% and .673893 at 50%. These are OOF results, not a new LB result;
no submission was made.

## Verification

The recorded run snapshots and checks all 40 primary pair-model memberships,
all ten pair prediction tables, and all five encoder training memberships.
The nested primary features are independently reconstructed from pair
predictions with absolute tolerance1e-12. Every training recipient excludes its
own nonanchor fold and the outer fold; every outer recipient uses only pair
models that exclude the outer fold. The encoder outer holdout sets exactly
match the fixed parent assignments and have no training/inner-selection
compound overlap. All 6,384 OOF keys, labels, and confidence bounds match both
old component and final ensemble. Input feature order, float32 conversion,
model parameters, and random seeds reproduce the original tuned readout
recipe, with only thread count reduced to three.

Artifacts include all 20 fitted models, `external-oof.parquet`,
`aligned-oof.parquet`, `comparison.csv`, `macro.csv`, and `summary.json`.
