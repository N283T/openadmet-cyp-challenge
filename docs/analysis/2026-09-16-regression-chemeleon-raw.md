# Regression CheMeleon raw2048 readout and ensemble placement

## Decision

Retain CheMeleon raw2048 + nested log2fc4 TabICL as a CYP2D6-specific full
ensemble-member candidate. Do not add it uniformly to all four endpoint
ensembles, and do not deploy both raw and PCA views everywhere. No test
inference, candidate packaging, or submission was performed.

Normal ensemble development now treats a new model as a complete equal-weight
member and reruns component selection. Incremental 10/25/50% additions are
reserved for deadline-defense diagnostics rather than routine exploration.

## Matched raw versus PCA readout

Run `20260916T060715739725Z-b0ddf3c8` uses the same adapted CheMeleon embeddings,
four nested out-of-fold log2fc predictions, TabICL 2.2 checkpoint, eight
estimators, outer UMAP folds, and train-only standardization as the PCA256 run
`20260912T071619253798Z-73a6fb7e`. The only representation change is retaining
all 2,048 standardized embedding coordinates instead of fitting fold-local
PCA256. This produces 2,052 final columns per fit and 6,384 keyed OOF rows.

Matched audit `20260916T061410497124Z-f652cf19` compares raw2048, PCA256, and
their fixed equal prediction mean without fitting weights.

| Endpoint | PCA256 | raw2048 | raw + PCA equal mean |
|---|---:|---:|---:|
| CYP1A2 | .744349 | .731283 | **.720223** |
| CYP2C9 | .501505 | .503371 | **.486256** |
| CYP2D6 | .903787 | **.866233** | .868608 |
| CYP3A4 | .451040 | .448192 | **.433319** |
| Macro | .650170 | .637270 | **.627101** |

The equal raw/PCA mean improves on PCA in all 20 endpoint-fold cells. Prediction
Pearson correlation ranges from `.8804` for CYP2D6 to `.9678` for CYP3A4;
residual correlation remains `.9253–.9613`. Raw therefore changes the readout
meaningfully, but does not constitute an independent error source.

## Placement as complete members

Authoritative peer audit `20260916T062147025088Z-1d1208a1` adds raw as one complete
equal-weight member of each preferred endpoint recipe. For CYP3A4 it uses the
reported cross-fold weighted top-two internal Caruana compression and never
duplicates a view already present in the current recipe. Initial audit
`20260916T061537607859Z-10752681` selected the distinct `nested_topk` row with
the same `k=2`; follow-up `20260916T061912834383Z-35d755e9` corrected that row
but double-counted the already-present PCA member for CYP2C9 in its dual-view
diagnostic. Both are superseded and are not used below.

| Endpoint | Current | + raw full peer | + raw and PCA full peers |
|---|---:|---:|---:|
| CYP1A2 | .687942 | .688116 | .690147 |
| CYP2C9 | .463103 | .463994 | .463994 |
| CYP2D6 | .828397 | **.826451** | .836187 |
| CYP3A4 | **.398316** | .402149 | .402149 |
| Macro | **.594440** | .595178 | .598119 |

Uniform adoption is rejected. Raw helps only CYP2D6 in this direct peer test;
including both CheMeleon views as separate peers is worse overall.

## Equal-member subset reselection

Run `20260916T061647720360Z-c8456edc` reruns the prior endpoint-wise exhaustive
equal-weight subset procedure after adding raw as the sixth top-level component.
For every endpoint and component count, the subset is chosen on the other four
outer folds and evaluated on the held-out fold. There is no fractional-addition
or continuous-weight search.

| Endpoint | Best count | Cross-fit ST-RAE | Full-OOF selected members | Raw selected? |
|---|---:|---:|---|---:|
| CYP1A2 | 4 | .682158 | Caruana, Jeremy frozen, Jeremy scores, Jeremy FT | No |
| CYP2C9 | 3 | .463103 | Jeremy scores, Jeremy FT, CheMeleon PCA | No |
| CYP2D6 | 4 | **.825571** | Caruana, Jeremy scores, Jeremy FT, CheMeleon raw | **Yes** |
| CYP3A4 | 5 | .403800 | Original five including CheMeleon PCA | No |

The CYP2D6 four-member composition is selected identically for all five held-out
folds. It improves the current two-member point estimate by `.002826`, although
the conditional paired bootstrap interval for current minus best spans zero
(`[-.005535, .011775]` when expressed as the current count-two delta versus the
best count-four result). Under the earlier smallest-non-inferior-count rule the
two-member recipe would remain selected. Given the present exploratory risk
preference, keep the stable raw-containing four-member recipe as a serious D6
candidate, not as a proven replacement.

A shared global subset does not support raw: cross-fit macro ST-RAE is `.596969`
at five members and worsens to `.597812` when all six, including raw, are used.
The raw benefit is endpoint-specific.

## Implementation

`scripts/external_tabicl_readout.py` now accepts `pca: 0`, fitting the
`StandardScaler` on outer-training rows and retaining every scaled embedding
coordinate. Saved preprocessing records explicitly contain `pca=None`, the
embedding output width, and the preprocessing mode. PCA configurations retain
their prior behavior.

## Follow-up: flattening the CYP3A4 ensemble

The scored CYP3A4 recipe still treats a weighted 52/48 combination of two
selected LightGBMs as one top-level Caruana member, then averages that
pseudo-member with four other components. This nested weighting is now legacy
behavior and conflicts with the preferred single-layer member interface.

Run `20260916T062945192364Z-884691d4` compares single-layer alternatives. The
two selected Caruana LightGBMs are each treated as independent members; all
members then receive equal top-level weight.

| CYP3A4 recipe | ST-RAE | Spearman |
|---|---:|---:|
| Current nested five | **.398316** | .825139 |
| Flat six: two selected LGBMs + frozen + scores + FT + PCA | .400327 | .820060 |
| Flat seven: flat six + raw | .402579 | .819940 |
| Flat four: remove both selected LGBMs | .409099 | **.827813** |
| Flat five: remove both selected LGBMs + raw | .410985 | .826376 |

Flattening costs `.002011` ST-RAE and improves three of five fold point
estimates. Removing the two LightGBMs entirely costs `.010782`; raw does not
recover that loss. Keep the flat six-member recipe as a structural baseline,
not a final ensemble choice. The currently scored vector remains unchanged
until a new package is explicitly prepared and submitted.

## Flat prediction bank and multiple ensemble methods

Run `20260916T064121913990Z-bad74b14` establishes the first regression flat-bank
interface. Its five eligible atomic OOF predictors are Jeremy frozen, Jeremy
scores+log2fc, Jeremy FT+log2fc, CheMeleon PCA+log2fc, and CheMeleon raw+log2fc.
The current endpoint recipe is a control only. Every fitted ensemble predicts an
excluded outer fold; Ridge, simplex regularization, and Caruana step count are
selected again within the four outer-training folds.

| Method | Macro ST-RAE |
|---|---:|
| Current legacy control | **.594440** |
| Interval-loss simplex | .598031 |
| Caruana | .598332 |
| Equal-subset selection | .598977 |
| Equal mean, all five | .601588 |
| MSE simplex | .602681 |
| Median, all five | .605050 |
| Ridge | .616500 |
| Positive linear stack | .620018 |

The initial bank does not yet contain the two separate CYP3A4 LightGBMs, which
explains part of the control advantage. This is a reusable interface and method
comparison, not support for retaining nested architecture.

Matched runs `20260916T064215152245Z-b62d33b1` and
`20260916T064220610853Z-33eaec8e` compare PCA-only and raw-for-PCA swap banks
with the dual-view bank.

Macro summaries are retained only as descriptive inventory. Bank composition
and ensemble method are selected independently for each target:

| Target | Best evaluated flat bank/method | ST-RAE | Current control |
|---|---|---:|---:|
| CYP1A2 | PCA bank / equal subset | .687942 | .687942 |
| CYP2C9 | PCA bank / Caruana | .463138 | **.463103** |
| CYP2D6 | raw-swap bank / Caruana | **.827973** | .828397 |
| CYP3A4 | PCA bank / Caruana | .407714 | **.398316** |

For CYP2D6, PCA-only interval simplex scores `.832718`; replacing PCA with raw
improves it to `.827988`, while including both gives `.828203`. Flat Caruana on
the raw-swap and dual banks scores `.827973`: raw receives positive weight in
all five outer fits and PCA receives zero when both are available. Equal-subset
selection uses scores+FT+raw in four of five outer fits. Raw is therefore a
supported D6 swap/addition. This does not imply a change for the other targets.

Run `20260916T065043779157Z-99517b89` then refits and exports the fixed
`class_quantile_0.65` and `morgan_quantile_0.65` CYP3A4 LightGBMs as separate
atomic OOF vectors. Runs `ef9f2b8d`, `f3c1c6ae`, `47d59d80`, and `d4bc47c7`
instantiate independent 1A2, 2C9, 2D6, and 3A4 bank configurations.

| Target | Target-bank best method | ST-RAE | Current control |
|---|---|---:|---:|
| CYP1A2 | Equal all | .690147 | **.687942** |
| CYP2C9 | Interval simplex | .463448 | **.463103** |
| CYP2D6 | Caruana | **.827973** | .828397 |
| CYP3A4 | Equal subset | .400035 | **.398316** |

The CYP3A4 equal subset is identical in all five outer selections: Morgan
quantile LightGBM, Jeremy frozen, Jeremy scores, Jeremy FT, and CheMeleon PCA.
It excludes the class-embedding quantile LightGBM and raw CheMeleon. The
remaining `.001719` gap to the legacy control is therefore measured after
atomic export, not caused by a missing bank member.

Future additions should append or swap atomic columns in the relevant target
bank and rerun this unchanged suite. Aggregate ensemble predictions remain
controls rather than recursively nested members.
