# Report-derived experiments: final aggregate result

Expanded global selection `20260912T094021204217Z-84c28266`; independent audit `20260912T094029976646Z-846181c4`.
The pool expands from 107 to 126 candidate identities. All four research directions and a combined-descriptor check completed; model outputs are preserved.

The comparison below uses leave-one-parent-fold-out weight fitting on the same retained OOF keys. Underlying models are reused: this remains a weight-stability check rather than newly nested model CV.

| Target | Previous global | Expanded global | Difference |
|---|---:|---:|---:|
| CYP1A2 | 0.674078 | 0.670794 | -0.003284 |
| CYP2C9 | 0.442189 | 0.443453 | +0.001264 |
| CYP2D6 | 0.825782 | 0.825436 | -0.000346 |
| CYP3A4 | 0.390980 | 0.391582 | +0.000602 |
| Macro | 0.583257 | 0.582816 | -0.000441 |

The macro movement is very small. CYP1A2 improves modestly; CYP2D6 barely moves. CYP2C9 and CYP3A4 regress slightly against the earlier global pool. No substantial focus-target breakthrough was found in this round.

Full-OOF weight-fitting scores are listed separately because the same labels select the weights:
- Previous: 0.573705
- Expanded: 0.572518

## New candidates selected in the full-OOF recipe

This is the fitted recipe, not evidence of an independent gain for each selected component. Complete weights, including five separate weight-crossfit recipes per endpoint, are saved in selected-models.csv.

### CYP1A2

| New candidate | Weight | Source |
|---|---:|---|
| `20260912T093227551001Z-ed3e4d03` | 22.00% | geometry_mean_std — Controlled ligand-only geometry additions to frozen600 plus nested log2fc |
| `20260912T093642111763Z-4ec73c9f` | 20.00% | joint_geometry_substrate — Matched simultaneous geometry substrate and Jazzy additions |
| `20260912T093601580048Z-ec2b7210` | 2.00% | independent_lgbm — Matched independent vs Hill observation MLP; all heads held-out; inner epoch selection and refit |

### CYP2C9

| New candidate | Weight | Source |
|---|---:|---|
| `20260912T092824781073Z-fef7788e` | 2.27% | report_tdi_control_newaux_v1 — Readout of completed fold-specific FT embeddings: report_tdi_control_newaux_v1 |

### CYP2D6

| New candidate | Weight | Source |
|---|---:|---|
| `20260912T093615055236Z-2a25a497` | 18.00% | report_tdi_arm_delta_v1_ft_lgbm — Fold-specific adapted external600 readout: report_tdi_arm_delta_v1_ft_lgbm |

### CYP3A4

| New candidate | Weight | Source |
|---|---:|---|
| `20260912T092909986920Z-c5d42521` | 2.00% | geometry_mean_std — Controlled ligand-only geometry additions to frozen600 plus nested log2fc |
| `20260912T093601580048Z-ec2b7210` | 2.00% | independent_direct — Matched independent vs Hill observation MLP; all heads held-out; inner epoch selection and refit |

## Branch conclusions

- TDI continuous potency and paired shifts produce small matched 1A2/2D6 LGBM gains but no macro gain in fixed additions.
- Hill-linked TDI response heads underperform matched independent heads in the frozen-feature MLP; the MLP control itself is weaker than existing models. This is a negative result for this implementation, not a refutation of the mechanism.
- Geometry features help CYP1A2 modestly, including the TabPFN comparison; CYP2D6 remains weak. Topology/basic-N and Gasteiger quantities are proxies, not pKa/microstate estimates.
- Real external substrate probabilities have little downstream benefit. This branch tests substrate supervision; atom-level SoM remains untested.
- Simultaneous geometry, substrate and Jazzy descriptors also fail to improve the prior manual ensemble at the prespecified fixed weights.

All fitted models, descriptors and OOFs are retained. No full-data test refit, submission or production-pointer replacement is included in this round.

Summary run: `20260912T094243016404Z-c8a4a010`.
