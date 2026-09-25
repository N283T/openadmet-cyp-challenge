# CYP2D6 regression: whole-member intermediate-depth screen

## Design

This screen follows the requested whole-model comparisons. The current D6
ensemble has two complete members: Jeremy score47 + nested log2FC4 TabPFN and
fold-safe Jeremy FT600 + nested log2FC4 LightGBM. Its equal mean is the control.
All candidates use the same parent-series five folds and 1,457 nonanchor D6
labels. No fractional weight sweep is used. Adding a model creates a third
equal member; replacing the FT member preserves two members. A multi-depth
block first averages its constituent readout predictions equally, then enters
the outer ensemble as one member.

CheMeleon depth-1 TabPFN and TabICL use fold-training-only standardization and
PCA256, with and without the four nested log2FC features. LightGBM readouts
for depths 1, 3 and 6 are from the matched parent-CV run `dcf9b1b8`. TabPFN
runs are `14683186` (embedding only) and `0fdef233` (+log2FC); TabICL runs are
`67f7eeec` and `d56c7302`.

The Jeremy depths come from the actual five fold-specific, label-safe FT600
encoders used by the current LightGBM member. Extraction run `227d97fa`
applies each encoder's frozen final atom projection to message-passing states
after steps 1, 2 and 4, then the original norm aggregation. Depth 4 matches
the saved final embeddings within `1.79e-7`. Run `39e028eb` refits the D6
LightGBM with the original fold-specific, inner-selected hyperparameters and
nested log2FC4. For the final-depth control it uses the exact saved embedding
files; its predictions reproduce the original member with zero difference.
This isolates the representation change without retuning the readout.

## OOF results

Lower ST-RAE is better. The atomic member screen is recorded in run
`49fea86a`.

| Complete readout or internal block | ST-RAE |
| --- | ---: |
| Current Jeremy score TabPFN | .860552 |
| Current Jeremy FT600 LightGBM (depth 4) | .853844 |
| CheMeleon depth 1 LightGBM | .899902 |
| CheMeleon depth 1 TabPFN / TabICL | .956701 / .985814 |
| CheMeleon depth 1 TabPFN / TabICL + log2FC4 | .917237 / .954958 |
| CheMeleon depths 1+6 / 1+3+6 internal mean | .892042 / .889916 |
| Jeremy FT depths 1 / 2 LightGBM + log2FC4 | .877675 / .870097 |
| Jeremy FT depths 1+4 / 1+2+4 internal mean | .852668 / .852017 |

| Outer ensemble recipe | ST-RAE | Change from current .828397 |
| --- | ---: | ---: |
| Current two members | .828397 | — |
| Add CheMeleon depth 1 LightGBM as third | .831715 | +.003318 |
| Add CheMeleon 1+6 block as third | .831091 | +.002694 |
| Add CheMeleon 1+3+6 block as third | .830326 | +.001929 |
| Add CheMeleon depth 1 TabPFN + log2FC4 as third | .842632 | +.014236 |
| Add CheMeleon depth 1 TabICL + log2FC4 as third | .848812 | +.020416 |
| Add Jeremy 1+2+4 block as third | .826553 | −.001844 |
| **Replace current FT member with Jeremy 1+2+4 block** | **.824937** | **−.003459** |
| Replace current FT member with Jeremy 1+4 block | .825637 | −.002760 |

The best observed replacement improves three of five folds, but fold 0 worsens
`.841136→.847793` and fold 3 is effectively tied. Conditional parent-group
bootstrap run `3b20fa56` gives the 1+2+4 replacement a paired ST-RAE delta
interval `[−.009955, +.002465]` across 2,000 replicates (1,185 groups),
with 87.75% of resamples better. The interval includes zero. It excludes
model-refit uncertainty and correction for inspecting several scenarios.

## Decision and limits

CheMeleon does not justify a third D6 member under the tested readouts.
Jeremy's intermediate states show a small possible benefit when **replacing**
the current FT member with the 1+2+4 block, but the evidence is insufficient
to change the submitted vector. The saved September 14 D6 column also has a
blind-leaderboard-informed affine placement to mean `3.57` and sample SD
`.90`; this OOF screen compares the underlying unplaced models. A promotion
would require full-training refits, test predictions, and a separate transfer
check of the replacement's ranking and placement. No candidate was packaged or
submitted.

The earlier [fractional prediction screen](2026-09-25-chemeleon-incremental-ensemble.md)
answers a different question and is not used for this decision. Compact OOF,
metrics and bootstrap results are tracked in
`artifacts/tracked/d6-depth-atomic-20260925/`.
