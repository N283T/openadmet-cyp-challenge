# D6 failure-linked comparison and prospective probe

The candidate tests independent D6 learning, inner-selected probability operating points and centered fold inference. It uses no historical score mixture and no prescribed positive count. A4 is copied from the last sent payload.

## Matched historical splits

| Model | Split | AP | Nested bootstrap MCC |
| --- | --- | ---: | ---: |
| historical_morgan_lr_c005 | umap | 0.2710 | 0.0843 |
| historical_morgan_rbf_c1 | umap | 0.2758 | 0.0777 |
| historical_morgan_complement_a1 | umap | 0.2898 | 0.0895 |
| fresh_d6_network | umap | 0.2939 | 0.1088 |
| historical_morgan_lr_c005 | scaffold | 0.3073 | 0.1086 |
| historical_morgan_rbf_c1 | scaffold | 0.3284 | 0.1246 |
| historical_morgan_complement_a1 | scaffold | 0.3303 | 0.0833 |
| fresh_d6_network | scaffold | 0.3188 | 0.1698 |

## Connection to prior failures

October 1 held the D6 call count at 49 while exchanging 21 calls each way. Live MCC fell .308635 to .104986. The called identities changed at fixed count.

The conditional multi-block AP strategy linked to the probe improved AP, but its historical-rate projected MCC did not improve. The AP-mass90 counterexample survived full-pool and scaffold checks, so neither correcting pool scope nor changing splits alone explains all failures. Full-pool and labelled-only curves are retained as sensitivity evidence; missing labels stay missing. Historical 6.533% is only a retrospective comparison.
The conditional strategy result is not an independent validation of the exact full-development-selected blind payload. The single-model strategy is shown separately.

## Inference sensitivity

| Model | Split | Full refit calls | Centered-fold calls | Changed |
| --- | --- | ---: | ---: | ---: |
| historical_morgan_lr_c005 | umap | 107 | 51 | 56 |
| historical_morgan_rbf_c1 | umap | 102 | 78 | 24 |
| historical_morgan_complement_a1 | umap | 133 | 178 | 45 |
| fresh_d6_network | umap | 158 | 77 | 81 |
| historical_morgan_lr_c005 | scaffold | 84 | 98 | 16 |
| historical_morgan_rbf_c1 | scaffold | 104 | 70 | 34 |
| historical_morgan_complement_a1 | scaffold | 67 | 244 | 177 |
| fresh_d6_network | scaffold | 215 | 173 | 50 |

## Reviewed payload

- Against last_sent, CYP2D6_is_TDI: 77 positives, 54 additions, 26 removals, 23 retained.
- Against last_sent, CYP3A4_is_TDI: 249 positives, 0 additions, 0 removals, 249 retained.
- Against best_D6, CYP2D6_is_TDI: 77 positives, 46 additions, 18 removals, 31 retained.
- Against best_D6, CYP3A4_is_TDI: 249 positives, 0 additions, 0 removals, 249 retained.

## Interpretation and limits

This is a submission candidate for a prospective experiment. It does not establish that the historic reversals are solved. The five-member blind aggregation is not an independent nested-CV evaluation of an identical bag; it is an inference choice whose sensitivity is explicitly reported. One Live observation tests the combined new pipeline and cannot isolate one cause. No individual hidden label is reconstructed. Uncertainty is documented for the submission decision rather than an automatic veto.
