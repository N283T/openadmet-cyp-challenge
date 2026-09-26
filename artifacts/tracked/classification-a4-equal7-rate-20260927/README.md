# Submitted A4 equal-seven full-pool rate audit

Compact scientific outputs from the following recorded classification runs:

| Run | Output |
|---|---|
| `20260926T150611714864Z-d10195fa` | Historical TabPFN member on all UMAP held-fold molecules |
| `20260926T150342241560Z-da6a7d98` | Four Jeremy members on all UMAP held-fold molecules |
| `20260926T150058222462Z-cf35fc2e` | Exact external Mordred/PXR tree member predictions |
| `20260926T150746748801Z-e74748cd` | Recomposition of the submitted equal-seven recipe; audit |
| `20260926T150803596078Z-4d566758` | Count sweep, group-bootstrap intervals and outer-fold choices |

The full-pool Parquet files retain compound IDs, endpoints, split and fold IDs,
and null labels. The equal-seven `score` ranks all held-fold molecules before
the supplied-label mask is applied. Saved labelled-only scores are reproduced
to Spearman `0.99999989` and mean absolute difference `0.0000147`; blind scores
and calls are verified exactly against the accepted package.

The 275-call candidate was recorded in `classification/20260926T150937578353Z-d4bf03a9`,
candidate `20260926T150945663844Z-90e7ad7b`, submission
`20260926T150948735129Z-6fc9d05b`, sent attempt
`20260926T151009796789Z-b6d309a1`.

See [the analysis](../../../docs/analysis/2026-09-27-classification-a4-equal7-rate.md).
