# Labelled-only Boltz TDI integration: ranking signal without a call-boundary gain

The user explicitly deferred complete-pool extraction because of effort and
storage. The authorized question is whether retained predictions contain enough
incremental signal to justify considering that investment later. No new Boltz
extraction, base-model fitting, candidate, or submission was performed. Regression
remains at its adopted checkpoint.

Diagnostic run: `classification/20261001T133225940252Z-502d403a`.
Report run: `classification/20261001T133723258661Z-7983eb18`.
[Complete report and retained tables](../../artifacts/tracked/cyp-boltz-tdi-diagnostics-20261001T133225940252Z-502d403a/report.md).
The fixed configuration is
[`tdi-partial-pool-diagnostics-v1.json`](../../configs/boltz-trunk/tdi-partial-pool-diagnostics-v1.json).

## Scope and checks

D6 has 1,497 labelled molecules; A4 has 3,584. These are the entire currently
covered training cohorts for their target-specific trunks, not the 6,146-molecule
common train/auxiliary pool. Respectively 4,649 and 2,562 molecules lack trunk
features. All 750 blind compounds already have features for each target.

The exact labelled incumbent composition is reconstructed from saved atoms.
Compound IDs, endpoint labels and UMAP fold assignments match the recorded
measurement and assignment inputs. No missing labels become negatives.
All 500 previously fitted readouts are reused; molecular-only controls remain
ancillary and are excluded from trunk selectors.

Single additions/replacements are selected using the other four saved OOF folds
by AP or target-prior MCC. A joint selector chooses the single placement and
diagnostic rate. Multi-model blocks use inner held-fold selection among those
four folds for block sizes 1/2/3/5, whole-block/individual-member additions and
10/25/50% rank mixing. Greedy multiple replacements preserve existing member
positions: D6 replaces only its PCA/raw block; A4 permits up to three distinct
replacements. Keeping the incumbent is always available.

Multiple trunk models here mean different pooling/readout recipes of the same
frozen receptor/ligand extraction, not new conformations or independently
computed trunks. All selection is conditional on saved OOF: it is not fully
nested base-model refitting. Other-fold base predictions may have been fitted
using the reported fold's labels.

## Primary result

MCC below uses the historical diagnostic fraction applied to labelled molecules
only and the historical working scored-population prior. It is not an estimate
from complete-pool call selection or a deployment rate recommendation.

| Target | Procedure | AP | AUROC | Diagnostic MCC |
|---|---|---:|---:|---:|
| D6 | Incumbent | .296860 | .580764 | .110109 |
| D6 | AP-selected single placement | .306124 | .591360 | .099878 |
| D6 | AP-selected block/mixing procedure | .318808 | .611664 | .104965 |
| D6 | MCC-selected single placement | .300191 | .590498 | .099878 |
| D6 | Joint single placement/rate | .299980 | .587428 | .090250 |
| A4 | Incumbent | .522733 | .805784 | .465398 |
| A4 | AP-selected single placement | .531759 | .810283 | .463866 |
| A4 | AP-selected multiple replacements | .534418 | .809185 | .431656 |
| A4 | MCC-selected single placement | .531162 | .809688 | .460802 |
| A4 | MCC-selected multiple replacements | .531621 | .810365 | .462335 |

None of the selected trunk procedures improves the incumbent's historical-rate
MCC point estimate. The MCC-selected no-trunk composition controls retain the
same pooled MCC as the incumbents on both targets. At other fixed diagnostic
rates, some trunk procedures improve relative to the same-rate control, but
rate/model selection does not produce an overall improvement: D6's joint
selector chooses 4% on four folds and 6.53% on one, and worsens held-fold MCC;
A4 chooses 33.2% on every fold and still worsens.

D6's AP block procedure actually selects a single LightGBM readout in four
folds: core pocket, Fe, and remaining HEM are pooled separately. It mixes that
rank with the incumbent rank at 50%. One fold selects a five-readout block.
Thus the larger AP gain is not evidence that five trunks are generally better.
Its 45 diagnostic additions include 17 positives, while 45 removals include 18
positives. The AP gain is .021949, but its conditional paired group-bootstrap
95% interval is [-.008438, .048718]; the MCC delta is -.005144 with interval
[-.061440, .042826].

A4's AP-selected multiple replacements add 153 diagnostic calls and remove 153;
26 known positives enter and 48 leave. The MCC delta is -.033742 with a
conditional interval [-.055335, -.010263]. Its AP increase therefore provides
poor evidence for the operating boundary. A4's single AP placement remains a
smaller, uncertain signal: MCC delta -.001532, interval [-.015694, .012801].

## Complementarity and interpretation

D6's core-pocket/HEM-mixed LightGBM, selected in the earlier single AP screen,
has prediction Spearman .423787 against incumbent ranks; rank-residual Pearson
is .782860 and squared-rank-loss correlation .376014. A4's pocket/HEM-separate
TabPFN with CheMeleon2048 and Mordred1613 has corresponding values .812896,
.898042 and .773043. The D6 readout provides a more distinct ranking, while A4
mostly overlaps the strong existing ensemble. Neither establishes useful call
swaps at the tested historical boundary.

Residuals are percentile ranks minus binary labels, not calibrated probability
errors. Their shared label term can inflate correlation. Raw Brier scores are
retained separately; all atom comparisons, per-fold correlations and trunk-pair
correlations are available in the compact tables.

The 2,000 paired UMAP-group bootstrap draws hold selected scores and calls fixed.
They do not rerank, repeat recipe selection, refit base models or adjust for
multiple comparisons. Some single-placement AP intervals exclude zero, but
that conditional ranking evidence does not establish prospective MCC utility.

## Decision

Preserve the D6 Fe-separated ranking signal and A4 single-member replacement
signal as exploratory findings. Do not prioritize completing 7,211 additional
target/compound complexes on these results alone: the tested call-boundary
procedures do not improve either actual incumbent, and multiple A4 replacements
are harmful in this diagnostic.

Keep full-pool extraction deferred under the user's stated preference. These
findings do not select a blind positive-call vector or authorize an upload.
The accepted 249-call A4 performance control and last-sent 275-call probe remain
distinct; their prospective rate failure is not repaired by this analysis.

`pixi run check` passes; the existing test suite reports 73 passed and 8 skipped.
All diagnostic outputs and informative negative results are retained with hashes.
