# Intermediate Chemprop depth: closed direction

## Decision

Close this line of ensemble work. Keep the submitted regression and
classification predictions unchanged. The user chose not to pursue further
intermediate-depth experiments after reviewing the complete-member screens.
No new candidate or submission was made from this branch.

## Evidence across endpoints

Intermediate CheMeleon states sometimes improve standalone readouts and
produce meaningfully different predictions. The relevant decision is whether
they improve the currently submitted ensembles as complete members or
replacements. Those tests did not establish a gain:

| Track / endpoint | Current OOF | Strongest observed intermediate-depth result | Assessment |
| --- | ---: | ---: | --- |
| Regression 1A2 ST-RAE | .670833 | .670757 replacing CheMeleon PCA TabICL with a 1+3+6 LightGBM block | Paired series-group interval crosses zero widely |
| Regression 2C9 ST-RAE | .454711 | .454711 when the equal-subset bank is reselected with depth candidates | New members never selected |
| Regression 2D6 ST-RAE | .828397 | .824937 replacing Jeremy FT with a 1+2+4 block | Paired series-group interval `[-.009955,+.002465]` crosses zero |
| Regression 3A4 ST-RAE | .394567 | .394567 when Caruana is reselected with CheMeleon depth candidates | New CheMeleon members never selected |
| Classification 2D6 target-prior MCC | .110109 | .115312 after adding depth-1 CheMeleon TabPFN to the selected core | AP/AUROC fall; labelled-only group interval `[-.016270,+.027342]` crosses zero |
| Classification 3A4 target-prior MCC | .465398 | .465398 with the 1+6 block as eighth member | Tied at labelled operating point; no supported gain |

The classification screen covers only labelled UMAP OOF, not the full held-out
molecular pool needed to select a submission call rate. Prior D6 labelled-only
OOF improvements have reversed on the blind leaderboard. The slight D6
classification point gain therefore does not justify a deployment change.

## Record disposition

Retain the feature extraction, readout code, run records, keyed OOF, and
negative results as reproducible research evidence. The detailed analyses
are [standalone CheMeleon depths](2026-09-25-chemeleon-intermediate-depth.md),
[regression D6 whole members](2026-09-25-d6-depth-atomic-ensemble.md),
[other regression targets](2026-09-25-other-target-depth-members.md), and
[classification whole members](2026-09-25-classification-depth-whole-members.md).
The earlier [fixed fractional ensemble screen](2026-09-25-chemeleon-incremental-ensemble.md)
is retained as historical exploratory evidence; its blend recipes are not
adopted. Existing candidates, submissions, and inference packages remain as
they were before this branch.
