# Why D6 classification is weak: structural coverage and label routes

Recorded diagnostics: `classification/20260912T125901882388Z-4d006374` and matched-pool follow-up `classification/20260912T130049427499Z-7f106e62`. All identities, measurements, UMAP assignments, observed auxiliary tables and selected ensemble OOF are declared frozen inputs. No models are trained and no official labels are rederived.

## The labelled populations are different

| Official endpoint | Labelled | Positive | Both direct/TDI pIC50 available | TDI pIC50 only | Direct pIC50 only |
|---|---:|---:|---:|---:|---:|
| D6 | 1,497 | 324 (21.64%) | 1,493 | 4 | 0 |
| 3A4 | 3,584 | 764 (21.32%) | 2,334 | 1,249 | 1 |

All1,250 3A4 calls without a complete pair are supplied negatives. D6 has only4 such negatives. Similar overall prevalence therefore conceals different label routes and different chemistry. Missing canonical direct potency means no usable value in this saved table; it does **not** establish that the assay was never performed. The raw-call provenance investigation determines the underlying rule. The4,649 D6 and2,562 3A4 training compounds with no supplied class call also have neither potency value in this prepared observed table.

Only259 compounds have both endpoint classification labels. Their positive rates are19.69% D6 versus40.93% 3A4, very different from the complete3A4 labelled population. A paired comparison must retain this prevalence difference.

## Structural coverage is poor in original UMAP OOF

Morgan radius2,2048-bit fingerprints are compared with exact Tanimoto similarity. Each labelled OOF molecule uses only labelled molecules outside its own UMAP fold as neighbor candidates. Test uses available training labels; no test label is read.

| Nearest labelled training neighbor | D6 | 3A4 |
|---|---:|---:|
| OOF median Tanimoto | 0.2985 | 0.3415 |
| OOF fraction Tanimoto≥0.5 | 1.20% (18 molecules) | 4.27% (153) |
| OOF fraction Tanimoto≥0.7 | 0% | 0.53% |
| Missing-class training query median, own fold excluded | 0.2958 | 0.3149 |
| Test median, all labelled training candidates | 0.4711 | 0.5397 |
| Test fraction≥0.5, all candidates | 47.47% | 64.67% |

The full-training test comparison has more candidates. The second run controls this exactly: for each fold, its test queries and held-out queries use the **same** labelled training pool. Even then D6 test median similarity ranges0.380–0.409 versus held-out0.292–0.317; test fraction≥0.5 ranges34.27–42.27% versus held-out0.35–1.83%. For3A4 the matched test median is0.479–0.515 versus held-out0.333–0.348, with47.33–55.60% versus2.90–4.89% at≥0.5.

Thus the UMAP structural extrapolation is substantially harsher than the test structural neighborhood, independently of candidate-pool size. This is structural evidence only; it does not predict an unseen test score. Positive and negative OOF molecules within each endpoint have almost identical median similarity, so poor D6 coverage is not uniquely concentrated in its positives.

## Neighbor labels do not prove D6 is uniquely noisy

We compare the nearest-label disagreement against500 permutations of supplied labels **within each UMAP fold and cohort**, preserving each fold's prevalence while holding neighbor identities fixed.

| Cohort | Endpoint | Observed disagreement | Prevalence-controlled null mean | Null95% range |
|---|---|---:|---:|---:|
| All endpoint labels | D6 | 31.26% | 33.90% | [31.70%,36.41%] |
| All endpoint labels | 3A4 | 31.42% | 33.59% | [31.85%,35.23%] |
| Same259 compounds / same neighbor graph | D6 | 32.05% | 31.61% | [26.25%,36.88%] |
| Same259 compounds / same neighbor graph | 3A4 | 44.79% | 48.29% | [42.08%,54.44%] |

Both complete populations have only a modest2–3 percentage-point agreement improvement over prevalence. D6 does not have worse raw disagreement than3A4. The shared cohort has extremely sparse coverage: median similarity0.2644 and no held-out nearest neighbor≥0.5. Its modest differences fall inside the permutation ranges. These diagnostics cannot identify measurement noise separately from activity cliffs or weak fingerprint relevance, and do not support calling D6 intrinsically noisier.

## Errors by observed assay regime

Predictions are the selected saved ensembles from `b0a05986`: D6 `all_plain_pfn_graph_mean` and3A4 `lgbm_balanced_mean`, with their original inner-derived thresholds. These are retrospective subgroups; the observed potency values are not inference features.

| Evaluated subset | n | AUROC | MCC | Positive recall | Negative FPR |
|---|---:|---:|---:|---:|---:|
| D6, both potencies available | 1,493 | 0.5665 | 0.0959 | 10.19% | 4.70% |
| D6, TDI potency only | 4 | Undefined | 0 | No positives | 0% |
| 3A4, both potencies available | 2,334 | 0.6950 | 0.2547 | 69.24% | 42.10% |
| 3A4, TDI potency only | 1,249 | Undefined | 0 | No positives | 12.97% |
| D6, shared259 compounds | 259 | 0.6131 | 0.1701 | See raw confusion | See raw confusion |
| 3A4, shared259 compounds | 259 | 0.6430 | 0.2573 | See raw confusion | See raw confusion |

Removing the additional3A4 negative route substantially lowers its discrimination relative to the overall selected ensemble. The gap also narrows on the same259 compounds. It remains a real gap on the both-measured populations, so cohort composition is part of the explanation, not a complete explanation. Independent restricted-cohort size-control training is handled separately.

D6's selected calls miss291 of324 positives in the both-measured population. Of293 molecules with measured `TDI pIC50 − direct pIC50` in[0.3,1),274 are supplied positives, but only24 are called positive: recall8.76%. Even the359 molecules with strictly positive conservative shift bounds (TDI lower minus direct upper>0) include276 positives of which244 are missed. This localizes many errors to a measurable-shift population; it is not only indecision around overlapping potency intervals. Its low overall AUROC also means changing one threshold alone cannot explain all of the failure.

Observed shifts≥1 are not uniformly supplied positives (D6 fifty of89;3A4 forty-two of103). Shift bins and conservative interval differences therefore remain descriptive regimes, never replacement label rules. Source-call quality/availability rules must be consulted before interpreting apparent exceptions.

## Artifacts and limits

`neighbors.parquet` preserves molecule and neighbor IDs, folds, similarities and supplied labels. CSVs preserve coverage,500-permutation null summaries, label-route counts and selected-ensemble errors. The follow-up saves per-fold matched test/OOF similarities and exact TP/FN/FP/TN counts by route. The five test fold measurements repeat the same750 compounds and are not independent observations.

The work identifies a substantial evaluation-distribution difference and an important3A4 negative-label route absent fromD6. It does not establish label errors, a mechanistic cause, or expected leaderboard improvement. No observed assay feature was added to a deployable model, and no submission was made.
