# TDI structural-alert augmentation

## Conclusion

Explicit TDI bioactivation alerts do not improve the current LightGBM baselines.
Several motifs are strongly enriched on their own, especially
methylenedioxyphenyl on CYP2D6, but Morgan fingerprints already encode their local
structure. Adding 12 alert counts plus `any` and total-count features mostly changes
threshold behavior without adding stable discrimination. Do not adopt this panel in
the classification candidate.

## Design

Run `20260913T093248333738Z-60af2657` refit four balanced LightGBM views on the
fixed five UMAP folds:

- Morgan radius 2 / 2,048 bits, with and without alerts.
- Morgan + Mordred, with and without alerts.
- Identical model parameters and fold-specific inner-MCC threshold procedure.
- Both challenge endpoints, including all 1,497 CYP2D6 and 3,584 CYP3A4 labels.

The fixed panel covers terminal alkyne, terminal alkene, furan, thiophene,
methylenedioxyphenyl, primary/secondary amine, hydrazine,
aniline/aminophenol, thioamide, thiourea, cyclopropylamine, and an
N-dealkylatable tertiary-amine proxy. These are transparent RDKit SMARTS
approximations of the categories in the archived OpenADMET article. They are not
claimed to reproduce its unavailable exact SMARTS implementation.

## Results

Pooled OOF results are below. Parentheses show the alert-minus-baseline change.

| View | Endpoint | fixed-0.5 MCC | inner-MCC | AUROC | AP |
|---|---|---:|---:|---:|---:|
| Morgan | CYP2D6 | .0860 (-.0075) | -.0175 (-.0928) | .5596 (+.0030) | .2666 (-.0016) |
| Morgan | CYP3A4 | .1791 (+.0000) | .1817 (+.0019) | .6629 (+.0044) | .3254 (+.0067) |
| Morgan + Mordred | CYP2D6 | .0861 (+.0082) | .0118 (-.0214) | .5309 (-.0024) | .2549 (-.0062) |
| Morgan + Mordred | CYP3A4 | .2881 (-.0177) | .2900 (-.0293) | .7639 (-.0011) | .4453 (-.0012) |

The paired-bootstrap analysis is run
`20260913T093847481702Z-4f05ab66` (2,000 within-endpoint compound resamples).
For the decision-relevant inner-threshold policy, macro MCC changes are:

- Morgan: -.0454, 95% interval [-.0720, -.0193].
- Morgan + Mordred: -.0253, 95% interval [-.0451, -.0069].

Thus the observed degradation is not just a tiny pooled-point fluctuation. The
small Morgan AUROC gains do not translate to a usable MCC improvement because the
inner-fold thresholds transfer poorly, most severely on CYP2D6.

## Univariate signal and assay-label caveat

The feature audit confirms real motif-level signal. On CYP2D6,
methylenedioxyphenyl fires on 21 compounds, 12 positive, for 2.64-fold enrichment;
furan is 1.87-fold and thiophene 1.41-fold enriched. On CYP3A4, the implemented
cyclopropylamine proxy is 2.52-fold enriched (29/54 positive).

The CYP3A4 challenge endpoint includes 1,249 preincubation-only compounds defaulted
to negative even though no IC50 shift can be measured. The source article excluded
those rows for its mechanistic alert analysis. We intentionally retain them here
because the model is evaluated against the challenge label. This lowers the full
CYP3A4 base rate to 21.3% and dilutes some article-level enrichments; it is another
reason a chemically valid alert need not improve challenge MCC.

## Decision

Discard direct alert concatenation for the current candidate. Retain the tested
feature generator and OOF/test predictions for possible future use as a weak,
mechanistically interpretable diversity component, but do not tune or select
individual alerts on these same outer-fold labels.
