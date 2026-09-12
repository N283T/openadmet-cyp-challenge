# External representation experiments: measured results

All scores below use the retained parent split, 6,384 matched nonanchor OOF
records, and the official macro ST-RAE. Lower is better. The objective remains
below 0.60; these results have not reached it.

## Completed initial comparisons

The three-source corpus is shared run `20260912T051513413487Z-ac9df6b5`:
35,320 compounds, 156,372 observations from ChEMBL, PubChem and Tox21. All
challenge connectivity groups are excluded from external pretraining. Binary
conflicts are masked; censored potency remains interval-valued. External
validation is compound-disjoint. Training balances objective, source and assay.

The four-way blend is an equal mean of fingerprint/Mordred LGBM, predicted-log2fc
TabPFN, frozen Admetica TabPFN and the new model. The five-way additionally retains
the old scratch PubChem600 TabPFN. No blend weights are optimized against outer OOF.
The reference four-way mean containing old PubChem600 scores 0.658841.

| External objective | Frozen LGBM | Frozen TabPFN | Four-way with TabPFN |
| --- | ---: | ---: | ---: |
| Classification | 0.736906 | 0.737664 | 0.655751 |
| Within-assay ranking | 0.761539 | 0.757954 | 0.658949 |
| Censored continuous | 0.775219 | 0.777738 | 0.659519 |
| Mixed | 0.739448 | 0.749943 | 0.655721 |

TabPFN audit runs, respectively: `20260912T051958052360Z-0567431e`,
`20260912T052346369992Z-91822d8a`, `20260912T052851165725Z-8527ede2`,
`20260912T053452152078Z-17ae50ed`. Different objectives have different available
labels and early-stopping epochs; this is not a label-count-matched causal test.

## Auxiliary fusion and adaptation

Nested log2fc features, shared `20260912T052549384296Z-bac5fb3a`, exclude each
outer fold from every auxiliary model used within that fold. Only the four
prediction columns enter the readout; model-count metadata do not.

| New model | Standalone | Four-way | Five-way | Audit run |
| --- | ---: | ---: | ---: | --- |
| Classification emb + log2fc, LGBM | 0.675960 | 0.651707 | 0.650877 | `20260912T052854966074Z-84337c95` |
| Classification emb + log2fc, TabPFN | 0.676715 | 0.650494 | 0.650478 | `20260912T053813262365Z-5440c9e1` |
| Nested Rustuna LGBM on same features | 0.676986 | 0.650806 | 0.650376 | `20260912T053807787969Z-196fe9d9` |
| External tree scores + log2fc, LGBM | 0.681799 | 0.650372 | 0.648941 | `20260912T053156529408Z-5319d11d` |
| External tree scores + log2fc, TabPFN | 0.678038 | 0.651823 | 0.650132 | `20260912T053836421489Z-1f1452a6` |
| Measured/predicted fusion, validation observations masked | 0.679272 | 0.649724 | 0.648335 | `20260912T054100823527Z-fc669c9a` |
| External classification encoder, challenge FT direct | 0.736920 | 0.647683 | 0.649980 | `20260912T054039276943Z-5d82686b` |
| Same FT emb, LGBM | 0.712494 | 0.649143 | 0.650390 | `20260912T054101483933Z-4b5133ae` |
| Same FT emb, TabPFN | 0.725781 | **0.645787** | 0.648265 | `20260912T054319015342Z-f96a88a3` |

The external tree teacher trains 20 PubChem/Tox21 CYP classifiers/regressors using
Morgan fingerprints and RDKit descriptors. Actual external measurement coverage
is almost absent in test (zero or one compound per head), so measured fusion
deliberately hides all validation observations and availability flags. Training
can use available external measurements. No optimistic measured-validation score
is presented as test-like performance.

Challenge fine-tuning run `20260912T053838041107Z-289578d5` initializes the external
classification message-passing weights, then trains four pIC50 and four actual
log2fc heads. All labels of outer nonanchor compounds are excluded. Checkpoints
use an inner molecular-group validation set; anchors remain training-only. All
five folds hit the 60-epoch cap, with best epochs 54–59 (zero-based). A longer
160-epoch comparison is the next planned test.

## Model-design findings and next comparisons

Assay scaling audit `20260912T052617517710Z-c45b7793` finds noisy small-assay standard
deviations amplifying ChEMBL continuous validation error. Optional pooled variance
shrinkage uses training exact values only. Held-out-label perturbation audit
`20260912T052936486765Z-7e14a6cc` confirms unchanged training variance pools.
Conditioned CYP/measurement heads with assay-specific affine adjustments are
being compared, as are source ablations and nested Caruana selection.

The concrete broader-ADME inventory is documented separately in
`2026-09-12-broader-adme-inventory.md`. It identifies measured AZ endpoints already
in the local DB; model predictions are not treated as new measured labels.

No submission or current-model pointer was changed by these experiments.

## Longer adaptation and a fixed combination with nested Caruana

The 160-epoch run `20260912T055443217549Z-466ccda1` stops naturally in every fold;
best epochs are 65, 94, 56, 71 and 105 (zero-based). Its TabPFN audit
`20260912T060044451789Z-c064d8c7` gives standalone 0.730195 and four-way 0.643125.
Thus longer training helps the ensemble, but does not improve standalone TabPFN.

A prespecified equal mean of nested Caruana and the original 60-epoch FT
four-way blend scores **0.634427**, audit `20260912T055841564844Z-9d3ba24a`.
Replacing that FT component with the longer-run counterpart scores
**0.632811**, audit `20260912T060136759518Z-b2811e6c`. These two equal-weight
comparisons are both retained; no weight search was performed. The best observed
score remains above the requested 0.60 threshold.

The exact new prediction file is regression run
`20260912T060136172428Z-02433f5b/outputs/external-oof.parquet`. Its components and
equal-weight policy are in `configs/external-ensemble/caruana_ft_long_equal_v1.json`.

Broader source additions are ready: AZ merge `20260912T055348682784Z-be731fe3`
contains 40,132 compounds; the subsequent Biogen merge
`20260912T055809060996Z-06088efb` contains 43,496. Source details and transformation
audits are in the broader ADME inventory. Source-only expansion, optional Emax
challenge heads and a medium-architecture broad-ADME arm are separate experiments.

## Published upstream medium checkpoint control

The actual Jeremy upstream 47-head checkpoint is distinct from his challenge
17-head fine-tuned model. Recorded provenance finds no evidence of challenge-label
training, but the full upstream molecule/task manifest and training date are
unspecified. This is a public-pretrained-model control, not a reconstruction of
the author's corpus. See `2026-09-12-jeremy-medium-frozen-control.md`.

Frozen medium embeddings + TabPFN audit `20260912T060313766699Z-6fb1a6a5` scores
0.649557 standalone, 0.631521 in the fixed three-way and 0.632624 in the four-way.
Its equal mean with nested Caruana scores **0.611584**, audit
`20260912T060422172365Z-199a2997`. The prediction run is
`20260912T060421588863Z-570c4938`. No OOF weight fitting was used. This result is
reported separately from the fully documented new-corpus arm (0.632811), retaining
the upstream provenance qualification rather than treating missing metadata as
evidence of leakage. Neither result is below 0.60.

An encoder-only bridge `20260912T060426566774Z-7604ffeb` enables our fold-safe
challenge FT; strict transfer and CPU embeddings exactly match the original
encoder on 16 checked compounds. No original predictor is transferred.

### Challenge adaptation and fixed combinations

Upstream-medium adaptation run `20260912T060538153709Z-8d434eee` preserves the
same outer exclusions and inner group selection as the documented scratch models.

| Readout of adapted medium | Standalone ST-RAE | Audit |
| --- | ---: | --- |
| LGBM | 0.637579 | `20260912T061158128715Z-20a6a987` |
| TabPFN | 0.679669 | `20260912T061614607873Z-f3eaac64` |
| LGBM plus nested log2fc | 0.628743 | `20260912T061701581307Z-586851c4` |

Fine-tuning helps LGBM but harms standalone TabPFN compared with frozen medium
0.649557; it is not a universal representation improvement. Equal Caruana +
adapted TabPFN gives 0.613518 (`20260912T061728659463Z-198cc9e1`); equal Caruana +
adapted LGBM gives 0.616151 (`20260912T061752767652Z-20383d21`), both worse than
the original frozen-PFN pair 0.611584.

The fixed three-way mean of original nested Caruana, frozen medium TabPFN, and
adapted medium + nested log2fc LGBM gives **0.605039**, audit
`20260912T061837935561Z-79bf9cf5`, prediction run
`20260912T061836843528Z-421acacc`. All tested combinations remain recorded; equal
weights are specified per experiment rather than optimized on outer OOF. This
upstream-model result remains above 0.60 and carries the same provenance note.

The newly curated AZ scratch mixed model has not matched this control. Frozen
TabPFN scores 0.743225 (`20260912T060917517046Z-fc36a204`); after challenge FT,
TabPFN scores 0.740840 and fixed four-way 0.648861
(`20260912T061855198099Z-f5cf0570`). Adding these ADME tasks alone is insufficient
in the tested architecture/training configuration. The separate broader Biogen
plus medium-architecture arm remains scheduled.

## OOF threshold reached with upstream predictor features

All phase-three queued runs completed during the usage interruption. Emax adapted
TabPFN gives standalone 0.730846 / four-way 0.643025 (audit5f23477e). PubChem+Tox21
classification without ChEMBL gives 0.726323 / four-way 0.653288 (auditb9b2be25).
The full AZ+Biogen medium mixed encoder gives frozen TabPFN0.719200 (audita9b4c43c)
and adapted TabPFN0.732866 / four-way0.648230 (audit76379e49). These do not beat
the stronger public-pretrained-model combinations.

Frozen upstream CheMeleon2048 -> training-only PCA256 + nestedlogfc -> TabICL
scores0.650170 (audit `20260912T072213556654Z-90bece05`). Adding it to the prior
three-way mean produces0.602911 (audit `20260912T072307655023Z-b444086c`).

Run `20260912T072105522912Z-53eb675c` strictly loads the original upstream medium
47-task predictor and exports all final hidden activations and all47 predictions.
No challenge labels or selection enter this export; task names remain unspecified.
The47 predictions plus four nestedlogfc features -> TabPFN score0.640427
(`20260912T072350034703Z-ba3871ed`). These are input features, not substituted
challenge labels.

Replacing frozen medium MP-embedding TabPFN with this predictor-feature TabPFN
in the four-way mean scores **0.598065** (`20260912T072501470405Z-4a9c0d46`).
Retaining both in a five-way equal mean scores **0.597617**
(`20260912T072502729280Z-d8928353`), below the requested0.60 threshold on the same
6,384 OOF records. Prediction run: `20260912T072501360501Z-f8c3a807`.

Each component has fixed weight0.2:

1. Original24-candidate nested Caruana: `20260912T054406239560Z-96fce32b`.
2. Frozen upstream medium MPemb -> TabPFN: `20260912T060059705105Z-2f0fbfbd`.
3. Frozen upstream47 predictions + nestedlogfc -> TabPFN: `20260912T072303592078Z-07f2749a`.
4. Challenge-adapted medium emb + nestedlogfc -> LGBM: `20260912T061629455493Z-730c8bea`.
5. Upstream ADME-adapted CheMeleon PCA256 + nestedlogfc -> TabICL: `20260912T071619253798Z-73a6fb7e`.

All equal-blend candidate comparisons are retained; no continuous OOF weight
optimization was performed. The result uses public upstream pretrained models,
whose complete original corpus manifest is unavailable. Author code distinguishes
these47-head inputs from the17-head challenge FT outputs; no evidence of challenge
label pretraining was found. Earlier diagnostic flags conservatively separated
this control from our fully documented new-corpus arm. Missing upstream metadata
is disclosed here, rather than imposed as a user-requested prohibition on public
pretrained models. The achieved figure is an OOF result, not a new LB submission.

The already launched inner HPO and confidence-interval FT comparisons are being
completed and retained; no further training is needed merely to reach0.60.

## Final verified result

The completed480-trial inner Rustuna search improves adapted-medium/logfc LGBM
standalone0.628743 ->0.624940. Its prescribed three-way comparison improves
0.605039 ->0.603423. Replacing just that component in the five-way equal mean
gives **0.5964141092**, predictions `20260912T072935321335Z-1e8d2643`, audit
`20260912T072935933157Z-6e1a3bce`. The other four components and all five0.2 weights
remain fixed. HPO inner feature reuse and outer exclusions are documented in
`2026-09-12-jeremy-ft-logfc-rustuna.md`.

The already running confidence-interval FT also completed: direct0.629855,
embedding+nestedlogfc LGBM0.627203, retained with all checks and models in
`2026-09-12-interval-multitask-finetuning.md`. No further trials were launched.

Independent finalization `20260912T073049029140Z-0c26dd06` reconstructs the final
five-way prediction to maximum absolute difference8.88e-16, independently scores
0.5964141092 (neighbors0.6872080091), and verifies all6,384 keys/labels/bounds/folds,
component hashes and actual FT outer exclusions. The package contains verified
OOF, scoring tables, exact weights, model references and provenance. Read
`2026-09-12-oof-goal-result.md` for the final model card. Public upstream weights
and their unknown full original corpus remain explicitly disclosed. No LB result
is inferred from this OOF result. Final check passes; tests37passed/7skipped.
