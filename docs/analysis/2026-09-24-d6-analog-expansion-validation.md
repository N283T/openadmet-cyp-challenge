# CYP2D6 TDI validation under the blind analog-expansion design

## Why the standard outer folds can miss the leaderboard

The [official challenge design](../../references/official/cyp/announcement/20260907T135802Z/snapshot.md)
selects the blind set by taking top hits from CYP1A2, CYP2C9 and CYP3A4 and
buying similar compounds. CYP2D6 was not used to select those anchors. All 750
blind compounds were then measured in dose response across four CYPs, while
training dose-response labels form a sparse matrix. Thus the blind D6 target
population is selected by other-enzyme chemistry and has much denser assay
coverage than the labelled D6 training population.

The recorded alias inventory `fe5cc874` finds 735 of 750 blind molecules in
the `OCNT-253` ID block and zero of 6,146 training molecules there. The ID
block is batch provenance, not a chemical feature or a class label. A separate
cross-fitted Morgan/descriptor train-versus-blind classifier had AUROC `.877`
(`4135649a`), consistent with a substantial distribution difference. That
classifier's density weights worsened AP-versus-control leaderboard concordance,
so they are not a transfer correction.

## Pseudo analog-expansion check

Run `1294f049` uses only training molecules to select 25 unique anchors per
1A2, 2C9 and 3A4 by measured direct-inhibition pIC50. It finds ten nearest
nonanchor training molecules per anchor by radius-2, 2,048-bit Morgan
Tanimoto similarity. The 750 anchor–analog pairs collapse to 608 unique
analogs. This approximation uses direct potency because the exact historical
hit-selection and in-stock purchase process is unavailable in the training
files. It does not use D6 labels for selection.

Only **164 of 608** selected analogs have a supplied D6 TDI label; 28 are
positive. Ranking the full outer-fold population before looking only at the
analogs is the wrong operating population for this question (the initial
`1294f049` decision table is retained but superseded). Run `46d1b9ba` instead
ranks the 608 pseudo-test molecules and calls the top 40, matching 49/750.
With a 7.06% scored-population prior, AP minus September 15 projected MCC is
`+.01275` on UMAP (95% grouped interval `[−.27619,+.32645]`) and `+.11497`
on scaffold (`[−.10108,+.31845]`). Only 3–9 of the selected calls per model
land on a supplied D6 label. The intervals are too wide to choose a model, and
the point estimates still do not predict the observed AP LB loss.

## Decision

This check establishes a concrete mismatch in validation geometry and label
coverage. It does **not** recover the leaderboard ranking. Reweighting or
filtering existing sparse D6 labels cannot create the dense shift measurements
made on the analog-expansion test batch. The September 15 submitted D6 vector
remains the only leaderboard-validated control. The four-swap 75/25 probe is
frozen as candidate `20260924T043703823086Z-51106137`, not uploaded; this
analog check supplies no positive reason to promote it.

A stronger validation would need either a genuinely dense, independently
measured D6 analog-expansion set or multiple carefully recorded leaderboard
probes. The historical board has seven accepted states from one hidden set;
it cannot by itself identify compound labels or tune a reliable selector.
