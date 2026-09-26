# CYP1A2 calibration probe on the current submitted bank

## Question and lineage

The current scored CYP1A2 column is the September 18 foundation-v8
ten-member equal bank carried into the September 24 and September 26 regression
submissions. Its Live ST-RAE is `.508530` and the full 750-row prediction
mean/sample SD is `4.412000/1.169409`. The [Interim target analysis](2026-09-26-regression-interim-targets.md)
finds a gap of `.0336` to tenth place. This experiment asks whether a small
additional location or spread change is supported without replacing the model.

Screen run `regression/20260926T121803762252Z-15616963` evaluates a fixed
5×5 grid: additive shifts `−.2, −.1, 0, +.1, +.2` pIC50 and positive scales
`.85, .925, 1, 1.075, 1.15` around each OOF cohort's prediction mean. The
**full_parent** cohort is the actual ten-member equal-bank parent OOF, 1,372
keyed labels. Full ten-member scaffold OOF is unavailable. The separate
**proxy_parent** and **proxy_scaffold** cohorts average three current bank
members (all-pool and dense-pool LightGBM, all-pool TabPFN), also on 1,372
labels each. Proxy results check the direction of a correction but cannot
validate the complete submitted ensemble on scaffold.

The initial run `regression/20260926T121540473710Z-a53fb089` incorrectly
named a parent-fold source snapshot as scaffold. Its apparent second-split
result is invalid as a scaffold comparison; it is preserved but superseded by
the corrected run with explicit split identifiers.

| OOF cohort | Identity ST-RAE | Scale .925, no shift | Change | Fold wins |
| --- | ---: | ---: | ---: | ---: |
| Full parent bank | .670833 | **.669712** | −.001121 | 3/5 |
| Three-member parent proxy | .678259 | **.675877** | −.002382 | 3/5 |
| Three-member scaffold proxy | .685148 | **.681406** | −.003742 | 4/5 |

Scale `.925` with no shift is the lowest pooled ST-RAE point in all three
screened cohorts. Shifts of magnitude `.1` or `.2` do not beat it. The
correction does not change within-cohort Spearman, but increases OOF MAE
slightly in all three cohorts: `.563786→.564865`, `.569381→.569797`, and
`.570106→.570574`. These are small, post-screen OOF improvements. The test
set was assembled with 1A2 active analogs, so the training-label activity
mix and optimal output spread may differ. No Live gain is inferred.

## Prepared files and decision

Packaging run `regression/20260926T121919456616Z-f6b2cea4` applies the
fixed, monotone rule `4.412 + .925 × (prediction − 4.412)` to the current
750-row scored 1A2 vector. Its mean stays `4.412`, sample SD falls
`1.169409→1.081703`, and all ranks are unchanged. Mean absolute change is
`.072628` and maximum change `.233792` pIC50. The run validates two files:
1A2-only against the current scored submission, and a combined file adding
the already prepared D6 `.88348→.90` calibration probe. CYP2C9 and CYP3A4
remain identical. The [compact outputs](../../artifacts/tracked/regression-1a2-calibration-20260926/README.md)
include both files and the screen metrics.

Combined candidate `20260926T121937881707Z-0fd32af2` was prepared as
submission `20260926T121941727761Z-5ec77801`, SHA-256
`41be58c13c6266cc3db322c4bc45175874025805b40f0d5518f59b395d538827`.
Its dry run validates the 750 identities, regression schema and private
submission profile. It has **not been uploaded**. This is an exploratory
calibration-only candidate: separate endpoint leaderboard results would allow
the 1A2 and D6 adjustments to be read independently, while the existing
September 26 submission remains the scored performance control. The OOF
evidence supports a limited local 1A2 sensitivity test, not an expected
leaderboard improvement claim.
