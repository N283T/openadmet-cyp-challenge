# CYP3A4 raw embedding replacement and simplification

## Decision

Use the fixed equal mean of seven top-level members as the first A4 development
candidate:

1. historical CheMeleon PCA256 + nested log2fc4 + Mordred TabPFN stack;
2. Jeremy600 PCA256 TabICL;
3. Jeremy600 raw600 TabICL;
4. Jeremy600 PCA256 TabPFN;
5. Jeremy600 raw600 TabPFN;
6. external-class + Mordred LightGBM; and
7. PXR/Admetica.

This removes fitted Caruana weights, keeps the model-addition interface explicit,
and was selected identically in all five held-fold-excluded searches. Its OOF
AUROC/AP/target-prior MCC is `.80578/.52275/.46540`, compared with
`.80687/.51391/.45467` for the current scored 80/20 vector. The grouped MCC
delta is `.01074 [-.00787,.02930]`, with 85.7% positive probability. It changes
18 of 249 blind-test calls. This is promising but not statistically secure; no
candidate or submission was created.

Keep a six-member fixed mean as the compact fallback. It removes Jeremy
PCA256 TabPFN from the full-data choice and scores `.80610/.52076/.45927`.
Four outer selections use this composition; one substitutes Jeremy PCA256
TabPFN for raw600 TabICL, so it is less selection-stable than the seven-member
choice.

## Matched raw600 readouts

Run `a459b098` fits outer and next-fold inner predictions for four raw600 views
on the original A4 UMAP folds, using the same eight-estimator TabICL/TabPFN
settings as their PCA256 controls.

| Representation/readout | PCA256 AUROC | Raw600 AUROC | PCA256 AP | Raw600 AP |
|---|---:|---:|---:|---:|
| Jeremy600 TabICL | .76892 | **.77520** | .46205 | **.49093** |
| Jeremy600 TabPFN | .76616 | **.77726** | .46464 | **.47998** |
| External-class TabPFN | .69437 | **.71688** | not previously reported | .37344 |
| External-mixed TabPFN | .68490 | **.69914** | not previously reported | .36127 |

Fold-normalized fixed-call evaluation in `d5adff22` confirms the Jeremy gains:
raw TabICL reaches `.77645/.49513/.38241` and raw TabPFN
`.77735/.48000/.38087` for AUROC/AP/target-prior MCC, versus
`.76814/.46188/.36698` and `.76738/.46599/.37315` for PCA256. Raw and PCA
predictions remain correlated (`.851` TabICL, `.856` TabPFN), but retain useful
boundary diversity.

Raw external-class, raw external-mixed and the earlier CheMeleon-parent raw
LightGBM are weaker standalone and are excluded by every size-seven selection.
Raw therefore helps the Jeremy representation but is not a universal PCA
replacement.

## Pool and simplification comparison

Run `d5adff22` first compares prespecified Caruana pools, then exhaustively
evaluates every fixed equal-weight subset of the ten-candidate raw/PCA bank.
Subset selection for each component count uses only the other four folds.

| Design | AUROC | AP | target-prior MCC |
|---|---:|---:|---:|
| Current scored 80/20 | .80687 | .51391 | .45467 |
| PCA5 AUROC-Caruana | **.81159** | .53078 | .45467 |
| Raw5 AUROC-Caruana | .80838 | .52446 | .44087 |
| Raw+PCA7 MCC-Caruana | .80823 | **.52723** | .45620 |
| Fixed equal mean, six selected | .80610 | .52076 | .45927 |
| Fixed equal mean, seven selected | .80578 | .52275 | **.46540** |

Replacing both Jeremy PCA views outright is harmful: `raw5` loses MCC despite
better raw standalone models. Keeping both raw and PCA lets the MCC-Caruana
pool improve every point metric slightly, while the fixed equal mean gives the
larger operating-point gain and simpler deployment rule. The ten-member
expanded Caruana never uses the three extra weak raw candidates beyond the
raw+PCA seven pool and produces identical predictions.

## Run disposition

- `5e5cac53`: failed before fitting due to doubly prefixed input names.
- `a459b098`: authoritative matched raw600 fit.
- `7eb1c3e5`: calculations completed but final mixed-dtype serialization failed.
- `4fdd549a`: complete intermediate, superseded by added member metrics.
- `d5adff22`: authoritative raw replacement and simplification audit.
