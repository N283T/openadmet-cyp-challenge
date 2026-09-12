# Focused frozen-feature readouts

Run `20260912T082257132357Z-bf656b4e` completes480 inner trials across two
targets, five folds and three feature/model views. It adds a small complementary
2D6 signal; it does not materially solve either weak endpoint.

The retained parent CV is unchanged. Per fold, a grouped20% inner holdout picks
one of16 Rustuna candidates; only the selected model refits on all outer-training
labels. Views are upstream47scores+log2fc4 → RBF-SVR, upstream600emb+log2fc4 →
RBF-SVR, and600emb+47scores+Morgan2048+log2fc4 → LGBM. Standardization for SVR
fits only inner training during selection and only outer training for final fit.
No Mordred columns are used despite the fingerprint input being a combined
Morgan/Mordred artifact. No test outcome or prediction is used.

| Model | 1A2 standalone | 1A2 25% addition | 2D6 standalone | 2D6 25% addition | 2D6 50% addition |
|---|---:|---:|---:|---:|---:|
| Current ensemble | 0.683335 | — | 0.834741 | — | — |
| Scores SVR | 0.726524 | 0.681159 | 0.846699 | 0.827849 | 0.827250 |
| Medium SVR | 0.730948 | 0.684975 | 0.879897 | 0.837009 | 0.845201 |
| Combined LGBM | 0.705509 | 0.681005 | 0.836353 | 0.827602 | 0.825598 |

Addition weights25% and50% were specified in the run code before results. The
complementary2D6 signal is present for both weights and two different learners.
For combinedLGBM50%, neighborST-RAE improves .832878→.805794, though Spearman
slightly declines .478540→.476267 overall. ScoresSVR25% improves Spearman to
.481066. These are still small effects, far from the leaderboard improvement
needed. 1A2's small mean improvement with25% combinedLGBM accompanies neighbor
worsening .726437→.740517; do not adopt it as a robust1A2 improvement.

Independent audit `20260912T082727776208Z-d5d17605` matches all2,829 keys, labels,
confidence bounds and neighbor flags per method, and recomputes scores with the
original official scorer. Conditional1,000-draw molecular-group bootstrap gives
2D6 combinedLGBM25% delta−.007139 (95% interval−.012019,−.002846) and50%
delta−.009143 (−.018541,−.000538). 1A2 intervals include zero. This is a bootstrap
of already fitted development predictions, not selection-adjusted uncertainty or
an untouched validation claim. The nested auxiliary inputs also do not rerun the
encoder within each inner readout selection; this optimism is documented.

Models, inner trial histories, selected parameter sets and memberships are saved
per endpoint/fold/view. `focus-oof.parquet` contains all three keyed methods;
`verified-comparison.csv` records exact metrics and conditional intervals.
The submitted ensemble remains unchanged; no submission was made.
