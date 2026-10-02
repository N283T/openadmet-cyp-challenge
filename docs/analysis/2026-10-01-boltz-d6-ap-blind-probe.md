# D6 AP-trunk probe: 21 blind swaps, inspection only

The user considered an exploratory D6 probe because previous OOF MCC designs
have not reliably predicted leaderboard directions, then explicitly requested
inspection without submission. This work freezes the concrete call vector for
review. It does not change the earlier diagnostic results or grant upload
authorization.

Selection run: `classification/20261001T135322663569Z-94d28c32`.
Package/review run: `classification/20261001T135736242152Z-113f153e`.
Candidate: `20261001T135841266638Z-eb390bc4`.
Prepared local copy: `20261001T135846483051Z-64386ab4`.
Payload SHA-256: `5bb177ecfbfabf4c57828fe4b63b181768209cd0b344ec4242e97c59c7e55621`.

[Complete changed-molecule review](../../artifacts/tracked/cyp-boltz-d6-ap-probe-20261001T135736242152Z-113f153e/report.md)
and [exact swaps with SMILES and neighbor diagnostics](../../artifacts/tracked/cyp-boltz-d6-ap-probe-20261001T135736242152Z-113f153e/blind-swaps.csv)
are retained with keyed scores, selection diagnostics, sensitivity tables and hashes.

## Exact recipe

Apply the previously diagnosed AP-block procedure to all development OOF, using
held-fold selection for block size, placement and rank-mixing strength, then
choose the final distinct readouts on all development labels. This selects one
fixed LightGBM classifier on protein/ligand single embeddings and separate
core-pocket, Fe and remaining HEM pair-embedding pooling blocks. The five raw
fold-model blind predictions are averaged before ranking across the 750 molecules.
Blend that rank 50:50 with the rank of the retained September 15 D6 ensemble.

Keep 49 positives to isolate score-recipe changes; this is not a new rate
optimum inferred from labelled-only OOF. The earlier AP .296860→.318808 is the
conditional held-fold procedure result, not a validation score for this exact
final blind vector. Full-development inner selection AP .324749 is descriptive
recipe-selection evidence, not independent confirmation.

## Blind changes

| Trunk rank weight | D6 positives | Retained | Added | Removed |
|---:|---:|---:|---:|---:|
| 0% | 49 | 49 | 0 | 0 |
| 10% | 49 | 42 | 7 | 7 |
| 25% | 49 | 34 | 15 | 15 |
| **50% selected** | **49** | **28** | **21** | **21** |
| 100% | 49 | 13 | 36 | 36 |

Only the selected 50% mixture is frozen as a candidate. Other weights are
call-change sensitivity diagnostics, not independently selected or validated
alternative payloads. The selected mixture exchanges 43% of the original D6
positive set. It is not confined to the old threshold: old ranks 1, 5, 7, 8 and
9 leave the positive set. Overall blend-versus-control Spearman is .870581;
trunk-versus-control Spearman is .524629.

Added compounds have mean similarity-weighted top-ten positive-neighbor fraction
.318208 versus .407410 for removals. However, added compounds have a closer
nearest positive on average (.496959 versus .415768), and a less similar nearest
negative (.345459 versus .379829). These conflicting summaries do not independently
validate the new calls or estimate hidden labels. The saved neighbor diagnostics
were reused without modifying their labels or definition.

Varying only the trunk fold inference while holding the historical control fixed
changes 18–28 calls relative to the five-fold-mean candidate. Added calls receive
3.71/5 positive votes on average; removed calls receive .71/5. This is a sensitivity
check, not five independent confirmation experiments or calibrated uncertainty.

## A4, validation and delivery state

A4 is copied exactly from last sent payload
`20260926T150948735129Z-6fc9d05b`: 275 positives, zero changed calls. It differs
by 26 calls from the earlier best-performing 249-positive control. This review
does not silently restore that earlier A4 vector. A later decision to restore it
would require a separate, newly reviewed payload.

Saved historical D6 scores reproduce all 49 last-sent calls; the last-sent and
performance-control D6 columns agree. All 750 ID/SMILES mappings, finite scores,
boolean predictions, frozen hashes and track validation pass. The local dry run
reports D6 42 changed calls, A4 zero changes, and `sent=false`. Its recorded
classification eligibility time is September 27, 2026 at 12:09:55 JST, already
past at inspection; no submission opportunity was consumed.

No upload, new training, new trunk extraction, regression change or automation
was performed. The initial selection attempt `06a3767e` stopped because a reused
ID frame omitted fold; the initial package attempt `04f926d3` stopped because
identity columns were compared in different orders. Both failures were noted,
retained and corrected in new runs; neither created a usable candidate.

The prepared ID remains for inspection only. Any future upload requires explicit
approval for its exact payload immediately beforehand.

## Subsequent composition: exact A4 best-control restoration

The user then favored restoring the best-scored A4 vector, considered a small
A4 trunk blend as an alternative, and favored retaining the reviewed D6 probe.
Prepare the exact A4 restoration as the recommended lower-risk companion; do
not add another A4 trunk experiment to this payload. Keep 49 D6 calls as a
comparison condition, not an asserted optimum. No upload was authorized.

New composition run: `classification/20261001T141004641753Z-76884efa`.
New candidate: `20261001T141033295310Z-d8c08df3`.
New prepared local ID: `20261001T141041310215Z-c7faa79f`.
New SHA-256: `67dbe95a6f9c307c5c8cea6746c1d993e08ecc854393052babf35be38bb6ba19`.

The D6 column is bitwise identical to the previously reviewed AP probe. The A4
column is bitwise identical to performance control
`20260924T050705396203Z-10e1d3ad`. Relative to the last sent payload, D6 changes
42 calls (21 added/21 removed) and A4 changes 26 (none added/26 removed).
Relative to the best performance control, A4 has zero changes. This composition
introduces no new model fitting, threshold search or inference.

Validation and local dry-run pass for all 750 rows; `sent=false`. Classification
cooldown is already elapsed. The earlier prepared payload remains retained and
unsent with its original A4 column. The new recommended payload also remains
unsent and needs explicit approval for its exact ID before upload.

[New frozen composition and provenance](../../artifacts/tracked/cyp-boltz-d6-ap-probe-20261001T141004641753Z-76884efa/report.md)
include target changes, both frozen IDs, dry-run and hashes.

## Authorized delivery; scoring pending

The user explicitly approved submitting the reviewed new payload. Prepared ID
`20261001T141041310215Z-c7faa79f` was sent once in attempt
`20261001T141324626858Z-cd9cca73`. A fresh dry-run reconfirmed the exact hash,
750 identities, D6 42 changes and A4 26 changes before the interactive send.
All three official pre-upload Live snapshots were saved. The API returned a
confirmed receipt; local attempt and submission status are `sent`.

Attempt creation was October 1, 2026 at 23:13:24 JST; the next local
classification eligibility time is October 2 at 11:13:24 JST. This uses one
classification submission opportunity. Regression was not submitted or changed.

One bounded post-receipt check at 14:14 UTC still shows the previous September
26 submission timestamp on MA, D6 and A4 Live boards. No score change is claimed;
receipt confirmation is distinct from completed scoring. No recurring monitor
was created. Compact delivery evidence is in
[`delivery.json`](../../artifacts/tracked/cyp-boltz-classification-delivery-20261001T141324626858Z-cd9cca73/delivery.json).
The earlier prepared probe with A4 275 calls remains unsent.

## Live result: reject the D6 blend

Timestamp-matched official snapshots at 14:28 UTC now show the October 1 14:13
UTC submission on all three Live boards. Recorded audit
`classification/20261001T143048820695Z-05fb009d` verifies the sent hash, all 750
identities, exact A4 best-control column, submitted timestamps and macro arithmetic.

| Endpoint | Previous MCC | New MCC | Delta | Previous rank | New rank |
|---|---:|---:|---:|---:|---:|
| D6 | .308635 | .104986 | -.203649 | 25 | 110 |
| A4 | .484306 | .529911 | +.045604 | 10 | 1 |
| Macro | .396471 | .317448 | -.079023 | 17 | 48 |

D6 precision falls .341206→.158841 and recall .383490→.192393. Thus both fall
despite keeping 49 total blind positive calls. The 21-call swaps fail on the
scored Live population. This is not a rate-only change. The earlier labelled-only
historical-rate MCC point estimate was negative, but did not resolve the severe
prospective magnitude. The AP gain did not establish useful Live call selection.

A4 exactly reproduces the prior best MCC, accuracy, precision and recall after
the bitwise 249-call restoration. Its recovery is confirmed separately from the
D6 failure. Retain the original D6 vector and restored A4 vector as performance
controls; reject this 50% D6 composition for current Live and preserve it as
scientifically informative negative evidence.

The upload was accepted/scored; model rejection is an experimental decision,
not an API rejection. No rollback upload, new training, extraction or regression
change was made. A future rollback requires explicit user approval for the
specific payload and cannot send before October 2 at 11:13:24 JST.

[Recorded Live metrics and provenance](../../artifacts/tracked/cyp-boltz-classification-live-20261001T143048820695Z-05fb009d/report.md)
retain the full comparison. Do not infer which individual hidden calls were
correct or any hidden-half generalization from these aggregate scores.
