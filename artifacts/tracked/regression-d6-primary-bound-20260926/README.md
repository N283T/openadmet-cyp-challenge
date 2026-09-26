# CYP2D6 one-sided primary-bound candidate

Compact evidence and deployment artifacts from the recorded runs listed in the
[analysis](../../../docs/analysis/2026-09-26-regression-d6-primary-bound-candidate.md).
The retained keyed OOF tables, group bootstrap samples, source snapshots and
input hashes remain in the immutable local runs.

| File | SHA-256 |
| --- | --- |
| `physchem-metrics.csv` | `5eee4aa09b4f86d4ac79214d15bb45ad69897649abe7170013b7884df74607f1` |
| `multitask-wide-metrics.csv` | `2ecc37d0447553003c268d5f359234f236bad54fb460591ca4d5906be21dbfa5` |
| `multitask-compact-metrics.csv` | `f54fa0e0bc3bb543a6a990970c2b9d5e146a6515ffc253f7ffa5ef2f5bef5226` |
| `upper-bound-metrics.csv` | `541cf568d78763683fe0d83d06ea7cba20284c2c57301b716ff728eba0f3d52c` |
| `blend-metrics.csv` | `888564b07ff8c08bb64bd5e9a75cd15c1a5732c97bc5021f83cb27830ef47330` |
| `fold-metrics.csv` | `7273998962d8334e520ae68ab62b666312b9d035cd9ee75186f56d2e729cf2a4` |
| `bootstrap-summary.csv` | `9d5228d6f402fc1483a275d3c1c41c1f14279efff50576bb70475eb9939667e7` |
| `fullfit-summary.json` | `77c1ec634f680ffd99ea3491a373e00cf14d2668f63b508261c63f14fd7620ea` |
| `d6-upper-model.txt` | `97ae004d1a1a792e73fed5373713e4725fc4df6f19c2e7adf3ea93a72aa66596` |
| `prepared-predictions.csv` | `adc7c5a5b600a39dd5d3f6353842e3a1ac41cfc00fc5e27aeb322e78888b08b4` |

The LightGBM booster predicts a centered residual; add the D6 direct-training
mean `4.784202026339585` to its raw output. Feature order is sorted `score_`
columns followed by sorted `emb_` columns, as implemented in
`scripts/regression_d6_primary_bound_fullfit.py`. The final payload additionally
matches the raw member's mean and SD to the September 24 submitted D6 column
and blends them equally. The exact source inputs and transformation are captured
in full-fit run `regression/20260926T042258195328Z-2027e07c`; replay run
`regression/20260926T042406809391Z-6b28df55` reproduces its candidate and
model hashes. Prepared submission `20260926T042337084384Z-c3578831` has not
been sent.
