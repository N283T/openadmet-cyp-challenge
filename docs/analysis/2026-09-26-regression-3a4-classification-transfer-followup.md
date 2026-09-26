# CYP3A4 regression: classification transfer follow-up

**2026-09-26 (JST).** The strong CYP3A4 TDI classification result motivated a target-specific regression model search. The regression response is direct-inhibition pIC50; classification predicts time-dependent inhibition. These are related but distinct endpoints.

## Scored control and prior transfer

The scored regression 3A4 column comes from the September 17 bank. Its 2,300-row, five-fold parent OOF ST-RAE is `.394567`; Live ST-RAE is `.431253`. The bank already includes a LightGBM arm trained on continuous TDI representations. That addition improved parent OOF `.398316→.394567` and Live `.444983→.431253` (RL0115). The classification probability was a weaker readout and received zero bank weight. Classification CYP3A4 TDI was strong on its own task (Interim MCC `.5023`, Live MCC `.5299`), but that is separate evidence.

## New experiments

| Approach | Parent OOF ST-RAE | Comparison |
| --- | ---: | --- |
| Scored 3A4 bank | .394567 | Control |
| Cross-fit residual correction using TDI OOF probability, half strength | .404798 | Worse; correction without probability also worse (.404153) |
| TDI embedding LightGBM arm alone | .432441 | Existing atomic arm |
| New TDI embedding CatBoost alone | .444844 | Worse than LightGBM arm |
| Raw 10% CatBoost + 90% scored bank | .395242 | Worse |
| Fold-moment-placed 10% CatBoost + 90% scored bank | .393689 | Improvement of .000878 |

The residual experiment used group-disjoint calibration and evaluation halves within each held regression fold, repeated five times. TDI scores were classification OOF predictions, with a different fold map. Its TDI coefficient changed sign across fits (range `−.0806…+.0985`), and primary ST-RAE worsened. Run `regression/20260926T132242507386Z-77f0d9d4`.

The new CatBoost used the same fold-specific 600-dimensional TDI representation and four nested log2FC features as the prior LightGBM arm. It was trained separately on each parent training fold, with no held-fold 3A4 labels in fitting. Raw standalone and fixed raw blends did not beat the scored bank on ST-RAE. Run `regression/20260926T132534292013Z-dc82a774`.

Label-free placement matched the new model's held-fold prediction mean and SD to the scored bank before fixed blending. This gave a small OOF improvement at 10% weight. It won four of five folds; fold 0 worsened by `.002319` while the others improved by `.001154–.002482`. Among 1,000 bootstrap replicates sampling 1,083 molecular-series groups, the 95% interval for the ST-RAE change was `−.001870…+.000269`; 94.6% of replicates improved. At 25%, the OOF gain was smaller and the interval wider. Placement run `regression/20260926T132707217977Z-5b4d52c8`; grouped audit run `regression/20260926T132934277827Z-030e6fff`.

## Decision

Keep the September 17 scored 3A4 model for the next submission. The useful classification transfer is already present in that bank. The additional CatBoost peer is directionally promising after placement, but its gain is smaller than the uncertainty of the group audit and has no blind validation. The placement rule would also need full-fit deployment before submission. Keep 3A4 unchanged in the prepared 1A2/D6 probe.

The first residual run `regression/20260926T132224927707Z-c4fd9af9` failed because the script used a nonexistent Polars DataFrame method. It produced no usable scientific result; the corrected run above is the evidence.
