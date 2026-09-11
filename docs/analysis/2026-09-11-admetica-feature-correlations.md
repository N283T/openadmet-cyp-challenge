# Admetica outputs as features: all-train and test correlations

Run `20260911T145924198895Z-373de784` evaluates the four unchanged Admetica
inhibition probabilities and four pre-sigmoid logits as candidate features.
This concerns model outputs, not their 300-dimensional hidden embeddings.

The population is all 6,146 non-test compounds, including those lacking pIC50,
and separately all 750 blinded-test compounds. Joins use compound IDs. Label
associations use all available pIC50 rows, including anchors. Test labels are
not available or used.

## Same-CYP association with predicted primary log2fc

| CYP | Train Pearson | Test Pearson | Train Spearman | Test Spearman |
|---|---:|---:|---:|---:|
| 1A2 | -.613737 | -.660438 | -.644507 | -.710589 |
| 2C9 | -.577133 | -.670166 | -.596823 | -.662586 |
| 2D6 | -.607585 | -.648804 | -.574202 | -.638959 |
| 3A4 | -.580699 | -.657269 | -.599159 | -.664555 |

Probability increases with inhibition, whereas log2fc decreases. The negative
sign is expected. On train, pred_log2fc uses the four eligible auxiliary models
excluding the recipient fold; test uses the ten-model mean. Thus these are
whole-population feature associations, but the train/test difference also
includes different averaging protocols and must not be attributed solely to
population shift.

## Association with pIC50 and other feature columns

Same-CYP probability vs all-labeled pIC50 Spearman:
1A2 .430014 (n1412), 2C9 .449820 (n1285), 2D6 .338083 (n1493),
3A4 .599930 (n2335). Cross-CYP feature/target associations are also retained.

All Mordred columns and all 2,048 frozen CheMeleon coordinates were compared
with each probability/logit, separately by population. Complete columns use
all molecules; incomplete descriptors use pairwise finite rows with at least
30 observations. Constants are omitted. Every row includes its sample size.

Examples of associations with fully observed Mordred columns (Spearman):
- 1A2 probability vs FCSP3 on train: -.514610; vs Mp on test: .651127.
- 2C9 probability vs ETA_beta_ns: train .446506, test .466101.
- 2D6 probability vs nBase on train: .365930; vs SaaCH on test: .451016.
- 3A4 probability vs FilterItLogS on train: -.441900;
  vs SlogP_VSA6 on test: .458804.

Top lists include some descriptors defined on small subsets (e.g. MINdsN for
125 training compounds); their correlations should not be compared as if all
6,146 compounds were used. Coordinate-wise embedding correlations depend on the
chosen representation basis. Neither low individual-coordinate correlation nor
these top lists prove a new feature cannot be reconstructed from the full
existing feature block.

## Outputs and interpretation

- `score-features.parquet`: keyed probability/logit/log2fc features, population,
  and fold; all 6,896 compounds.
- `score-correlations.csv`: all cross-CYP score/log2fc and probability pairs,
  Pearson and Spearman, train/test separately.
- `all-label-correlations.csv`: every Admetica score vs every labeled pIC50 endpoint.
- `all-feature-correlations.parquet`: exhaustive per-column feature correlations.
- `top-feature-correlations.csv`: top five absolute Spearman associations per
  population, score and feature family, with sample sizes.
- `score-distributions.csv`: mean, SD and quantiles for probabilities/logits.

Retain the four output scores as a compact feature-addition candidate. Probability
and logit encode the same underlying ordering, so eight columns are not eight
independent signals. A matched feature-addition experiment is required to measure
incremental utility. This run performs no downstream fitting or submission.
