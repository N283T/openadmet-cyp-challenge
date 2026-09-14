# Morgan fingerprint attribution in the retained TDI classifiers

## Conclusion

The two endpoints use Morgan information differently. CYP2D6 has weak, sparse
univariate evidence and benefits most from radius-2 environments. CYP3A4 has many
reproducible associations and is driven more strongly by radius-0/1 local
chemistry. Neither endpoint is well described by a short list of the largest
logistic coefficients, and hashed bits must not be interpreted as unique
substructures.

## Scope and design

Run `20260913T124832303603Z-67007a7a` regenerates Morgan radius-2 2,048-bit
features and analyzes both labeled endpoints on the fixed UMAP folds. It combines:

- Univariate label enrichment with Fisher tests and Benjamini-Hochberg correction.
- Fold-specific coefficients from the exact weak logistic recipe (`C=0.05`).
- Heldout feature contributions and gain from the saved Morgan LightGBM component.
- Nested radius and top-bit ablations, with feature selection repeated inside each
  outer training partition.

Run `20260913T125041600108Z-b880dd05` extracts non-ubiquitous bits whose enrichment
direction, logistic sign and LightGBM usage are stable across folds. This explains
the Morgan components, not the full heterogeneous live ensemble. In particular,
the live D6 rank is a 2:1 blend of the prior 11-component ensemble and the Morgan
logistic rank.

## CYP2D6

Only two non-ubiquitous bits survive FDR 5%. Bit 1654 occurs in 105 labeled
molecules but only four positives (odds ratio 0.133, q=0.00045). Bit 1872 occurs
in 55 molecules with one positive (odds ratio 0.064, q=0.0348). No positively
enriched bit survives multiple-testing correction. Logistic and Morgan-LightGBM
nevertheless share 33 of their 50 highest-priority signed/gain bits, showing some
model-level agreement despite weak individual tests.

Nested ablation points to extended environments rather than atom identity or
first-neighbor chemistry:

| Feature set | AUROC | AP | Projected target MCC at live rate |
|---|---:|---:|---:|
| All radius-2 Morgan bits | 0.5729 | 0.2710 | 0.0384 |
| Radius-2 occurrences only | 0.5895 | 0.2857 | 0.0611 |
| Radius-1 occurrences only | 0.5594 | 0.2520 | 0.0125 |
| Radius-0 occurrences only | 0.5239 | 0.2393 | 0.0253 |
| Top 50, selected within training fold | 0.5314 | 0.2416 | 0.0210 |
| All except training-selected top 50 | 0.5888 | 0.2777 | 0.0519 |

The counterintuitive top-bit result matters: large coefficients are not a compact
mechanistic signature. Correlated and fold-specific weak bits carry more useful
ranking information, while some largest coefficients overfit the small endpoint.
The radius-2 result is a post-analysis hypothesis and needs a prespecified rerun
or prospective LB test before model adoption.

## CYP3A4

CYP3A4 has 110 non-ubiquitous bits at FDR 5%, compared with two for D6, and 232
bits meet the cross-fold/cross-model consensus filter. Logistic and Morgan
LightGBM overlap on 31 of their top 50 bits. Representative strong associations
include bit 1114 on the positive side (28.1% positive versus 21.3% endpoint
prevalence; odds ratio 1.61, q=0.000086) and negative bits 1911, 294 and 667
(positive rates 5.8%, 12.0% and 10.4%; odds ratios 0.21, 0.43 and 0.39).

The radius pattern is the reverse of D6:

| Feature set | AUROC | AP | Projected target MCC at live rate |
|---|---:|---:|---:|
| All radius-2 Morgan bits | 0.6562 | 0.3104 | 0.1879 |
| Radius-0/1 occurrences | 0.6837 | 0.3428 | 0.2130 |
| Radius-1 occurrences only | 0.6723 | 0.3307 | 0.2099 |
| Radius-2 occurrences only | 0.6076 | 0.2719 | 0.1264 |

This supports a broad local-chemistry signal for A4 and suggests that extended
radius-2 environments add variance rather than useful endpoint specificity in
the weak linear readout.

## Interpretation

Morgan bits are hashed buckets. A bit can represent several atom environments,
and the recorded fragment strings are examples observed in this dataset rather
than a unique SMARTS definition. The frequent top bits also include generic atom
environments. Chemical claims therefore require collision-free sparse Morgan
identifiers or explicit fragment features before naming a motif.

The actionable result is representational rather than mechanistic: test a
prespecified radius-2-only D6 view and radius-0/1 A4 view, preferably with unhashed
or larger fingerprints as a collision control. Do not construct manual alerts
from the current top-bit table, and do not replace the live classifier from this
diagnostic alone.
