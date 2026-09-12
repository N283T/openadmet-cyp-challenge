# CYP2D6 focused round: aggregate result

Selection `20260912T102554362041Z-7f98ba4d`; audit `20260912T102604139797Z-e5518a36`.
Expanded pool: 174 candidate identities.

The comparison holds each parent fold out when fitting ensemble weights.
Base OOF models are reused, so this checks weight stability, not fully nested model CV.

| Endpoint | Previous | Expanded | Difference |
|---|---:|---:|---:|
| CYP1A2_pIC50_direct_inhibition | 0.670794 | 0.670794 | +0.000000 |
| CYP2C9_pIC50_direct_inhibition | 0.443453 | 0.443453 | +0.000000 |
| CYP2D6_pIC50_direct_inhibition | 0.825436 | 0.825952 | +0.000515 |
| CYP3A4_pIC50_direct_inhibition | 0.391582 | 0.391582 | +0.000000 |

Other endpoint predictions match the frozen prior pool within 1e-12.

## New components in full-OOF fitted CYP2D6 recipe

These weights use the full OOF labels and are not an independent validation result.

| Run | Method | Weight |
|---|---|---:|
| `20260912T101411999351Z-601b7cb4` | d6_specialist_d6_direct_fresh_scratchrates_v1_lgbm | 6.38% |
| `20260912T101606024008Z-70d8bdfa` | d6_specialist_d6_direct_fresh_v1_lgbm | 2.13% |
| `20260912T101831048657Z-478527aa` | d6_joint_p05_lgbm | 8.51% |
| `20260912T101831048657Z-478527aa` | d6_joint_p30_lgbm | 4.26% |
| `20260912T102254368746Z-570cbe07` | d6_joint_p00_lgbm_refit | 2.13% |
| `20260912T102337278051Z-f89acb5a` | d6_specialist_d6_aux_seed17_v1_lgbm | 4.26% |
| `20260912T102357053293Z-e9975cb3` | d6_scope_mean_lgbm | 2.13% |

Summary run `9ee864b2`. First summary attempt `d4b9eda7` required bitwise
identity for other endpoints; harmless floating-point differences triggered it.
The completed summary verifies equality within 1e-12. All four endpoint pools
include the full-scope specialist's direct outputs; other new models are D6-only.

The expanded full-OOF fitted macro is 0.572292, versus 0.572518 previously.
Weight-crossfit macro is 0.582945, versus 0.582816: no improvement. The original
pool is frozen, the unconverged kernel run is explicitly excluded, and all other
completed variants are retained without score-based filtering.

Separately, the prespecified 50% addition of the fixed four-scope LGBM average
scores D6 0.822232 versus 0.825436. Canonical 6,384-row development predictions
are preserved in run `20260912T102818791182Z-1aec6ce1`, with other endpoints
exactly inherited. Conditional 2,000-draw molecular-group bootstrap delta interval
is [-0.009087,+0.002486], so this small point improvement is unresolved. It is not
a selected production pointer or a new submission.
