# Classification label provenance and official baseline audit

## Findings

The retained canonical classification labels agree exactly with the supplied booleans. There is no evidence here of a label inversion, broken alias join, or conflicting-label merge that explains weak CYP2D6 performance. This does **not** establish that insufficient data causes the weak ranking. The current submission is below the official TabICL classification baseline on both endpoints; the D6 gap to TabICL is small compared with the gap to the leading submission.

Recorded raw audit: `classification/20260912T130054345728Z-e960da8e`.
Recorded potency-rule coverage: `classification/20260912T130537108031Z-1edf91e1`.

## Raw-to-canonical reconciliation

The audit follows `(source_file, molecule_name, raw_smiles)` through the saved canonical alias table, verifies canonical IDs equal SHA256 of standardized SMILES, and checks all retained supplied labels against canonical measurements. The official scored endpoints have 5,081 labels: D6 1,497 with 324 positives; 3A4 3,584 with 764 positives. There are zero value disagreements, zero conflicting labels, zero unresolved alias keys, and zero retained labeled groups combining different raw SMILES or names.

The 10,162 retained labeled source occurrences include repeated labels in Emax and TDI files and collapse to the 5,081 canonical labels. This count excludes **one labeled source occurrence**, not two: Emax-only `OCNT-2329004`, raw SMILES `N.N.[Cl][Pt][Cl]`, has D6=False and missing 3A4. Its two endpoint rows are preserved in `quarantined-label-rows.csv`, but only one has a label. Molecular standardization quarantined this compound. Raw Emax therefore has 1,498 D6 labels, while raw TDI and canonical data have 1,497. The first audit attempt `dc511012` failed an overly strong assumption that every raw row has a canonical ID; the completed audit explicitly accounts for this existing quarantine. No source data were changed.

## Supplied labels versus the published potency rule

The official challenge snapshot describes a greater-than-twofold potency shift, inferred positives with direct pIC50 below 4 and TDI pIC50 above 4.301, and assigned negatives when both potencies are below 4. The audit uses the exact log10(2) threshold and requires both potency values to reconstruct a label. All 1,088 supplied positives agree with this reconstruction; zero defined reconstructed labels disagree with supplied labels.

Some supplied negatives cannot be derived by this simplified public rule:

| Endpoint | Missing direct potency | Missing TDI potency | Both available, rule undefined | Total unresolved negatives |
|---|---:|---:|---:|---:|
| D6 | 4 | 0 | 18 | 22 |
| 3A4 | 1,249 | 1 | 159 | 1,409 |

These rows remain valid supplied classification labels; the audit does not replace or discard them. The source definition is insufficient to independently reconstruct every negative. In particular, 3A4's labeled cohort includes many compounds without measured direct potency, whereas D6 mostly has paired potency measurements. This establishes a difference in label/cohort coverage, not a demonstrated mechanism for model performance. It would be incorrect to claim the simple rule reproduces every supplied label.

## Live leaderboard evidence

Snapshots fetched at approximately 2026-09-12 12:56 UTC using `cyp leaderboard`: MA `738ff4c9`, D6 `e4ea4a7a`, 3A4 `64e5e66c`. Our current row is timestamped 2026-09-12 12:27 UTC. Values below are bootstrap mean MCC and reported bootstrap standard deviation, rounded for display.

| Endpoint | Current N283T MCC (std) | Official TabICL MCC (std) | Current rank / TabICL rank |
|---|---:|---:|---:|
| D6 | 0.176698 (0.084810) | 0.199337 (0.090631) | 43 / 35 |
| 3A4 | 0.375036 (0.044817) | 0.450597 (0.052491) | 45 / 9 |
| Macro | 0.275867 (0.046806) | 0.324967 (0.053291) | 45 / 21 |

D6 rank 1 is Asidsal11: **0.430129 MCC, std 0.102652**. Its precision/recall are approximately 0.778816/0.267011, versus our 0.238594/0.231328 and TabICL's 0.273773/0.233414. The leading submission demonstrates that substantially higher D6 LB performance is possible on the challenge. It does not identify which training data or method causes that advantage. Independently reported standard deviations cannot support a paired significance claim about either gap.

Our 3A4 precision/recall are 0.452306/0.821323 versus TabICL 0.682377/0.504553, indicating quite different operating points. This is relevant counterevidence against explaining the entire shortfall solely as insufficient learned ranking.

Historical snapshot evidence in `docs/submission-baseline.md` records the September 10 submission at D6 MCC 0.001644, 3A4 0.430653, macro 0.216149. The current submission improves D6 substantially but worsens 3A4. The old row is a historical snapshot, not a simultaneously visible live leaderboard row.

## Official scoring implementation

The archived public evaluation implementation is pinned to `OpenADMET/CYP-Challenge-Tutorial` revision `858ae63ce79934113bccdb7fc65467de5f7b1935`. The audit additionally archives its `evaluation/utils.py` (SHA256 `c6ddfc1e8d8cf26c77567125680da8b9656bfa43fca31f62f0e333e21aba7cb6`). The public configuration describes itself as ported from the backend. The current public Space frontend delegates scoring to AWS Lambda; the executing private backend revision is not independently exposed by this snapshot.

The public evaluator right-merges predictions onto ground truth by `Molecule_Name`, sorts names, and filters missing ground truth separately for each endpoint. A missing prediction for a labeled row is an error. Supplied classification labels and predictions are cast to bool; labels are not recomputed from potency at scoring time. MCC uses sklearn `matthews_corrcoef`, including its zero result for degenerate cases. Precision, recall and F1 use `zero_division=0`.

There are 1,000 bootstrap samples with replacement, each the size of the endpoint's labeled set. RNG seed is 0. Reported endpoint values are the mean of bootstrap metric values; standard deviations use pandas sample standard deviation. Macro metrics are calculated by averaging endpoints within each bootstrap replicate, then summarizing replicates. Macro standard deviation is not the mean of endpoint standard deviations. These standard deviations are neither standard errors nor 95% confidence intervals. A pooled local OOF MCC is therefore not numerically identical in definition to the reported bootstrap mean MCC, in addition to evaluating different compounds.

## Historical and external model provenance

The old submitted model's OOF aggregate is `243a1e23`, threshold selection `c95bbbfc`, inference `343eb5ec`, prepared submission `68969a86`, and accepted attempt `f2dbb957`. Actual fold readout configuration (for example `492cf5f6`) is `variant=stack_mordred`, not merely the generic `features=pretrained_embedding_only` string also retained in that configuration. `descriptor_readout.py` concatenates Mordred with 260 stack features. Stack `fbc2b2ea` explicitly specifies frozen CheMeleon PCA256 plus four nested predicted primary-screen log2FC features, shared UMAP `edba6fc0`, and auxiliary fits excluding the outer fold and each corresponding inner fold. The readout is TabPFNClassifier 2.6 with eight estimators. This is a code/config provenance inspection, not a fresh numerical reconstruction of every historical fitted model.

The old D6 threshold 0.135978 produced **686 of 750 test positives**. Its OOF AUROC was 0.542860; the corresponding 3A4 OOF AUROC was 0.793282. Thus the historical D6 failure included both weak OOF ranking and a severe submitted operating-point problem. It should not be reduced to data size alone.

The current external-only 600-dimensional encoders `a3ea2ec3` (classification) and `2c0a0953` (mixed objectives) use dataset `ac9df6b5`: source-separated ChEMBL potency/threshold labels, PubChem AID1851 activity, and Tox21 AC50/activity observations. Preparation preserves assay IDs, masks conflicting binary labels, and excludes challenge connectivity. No matched paired direct-versus-preincubation potency-shift target is constructed by these preparation scripts. Therefore these experiments demonstrate transfer from general CYP assay supervision, not already-tested transfer from an external dataset reproducing the challenge TDI-shift definition. This is a coverage limitation, not proof that better matched external labels will solve the problem.

## Interpretation

Evidence supports weak D6 ranking in the older model and different classification operating points across submissions. It rejects a broad retained-label corruption explanation. Neither the raw audit nor the small gap to official TabICL establishes a sample-size ceiling or an unavoidable property of D6. The much stronger top D6 result, incompletely reconstructed negative-label rules, distinct D6/3A4 cohort coverage, and absence of matched external TDI-shift supervision all remain relevant counterevidence to a single-cause explanation.
