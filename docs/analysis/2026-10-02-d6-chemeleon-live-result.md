# CheMeleon-primary candidate: small recovery versus last sent, below historical best

**Identity diagnostic:** [the exact model top 49](2026-10-02-d6-chemeleon-top49-overlap.md)
shares 22 with the historical best and replaces 27 each way. The actual 175-call
payload shares 43, adds 132 and removes six. No counterfactual Live score is known.

Live leaderboards downloaded October 2 at approximately 21:52 JST show the
participant row timestamp **2026-10-02 12:08 UTC** in MA, D6 and A4. This matches
the [authorized 175-call submission](2026-10-02-d6-chemeleon-submission-receipt.md),
prepared `20261002T092348406577Z-4ea606c1`, attempt
`20261002T120839637505Z-f9a3d830`. The submission is scored, not just received.
All before rows are timestamped 2026-10-01 14:13 UTC.

| Live metric | Immediately previous | Current | Change | Current rank |
| --- | ---: | ---: | ---: | ---: |
| D6 MCC | .104985814 | .117740060 | +.012754246 | 105 |
| A4 MCC | .529910537 | .529910537 | .000000000 | 1 |
| MA MCC | .317448175 | .323825299 | +.006377123 | 48 |

The historical best D6 MCC is .308635; this candidate remains substantially
below it. D6 precision falls .158841077 -> .124688496, while recall rises
.192392837 -> .421533325. Calling 175 instead of 49 positives recovers some
recall without establishing a worthwhile MCC gain over the historical best.
A4 retains its exact best-scored 249-call vector and unchanged metrics. MA
equals the arithmetic mean of the two endpoint MCCs.

The completed matched refit favored this candidate's trained decisions on
observed development outcomes. That advantage did not produce superiority to
the historical best on Live. The earlier narrow-rank disadvantage and absent
weak-primary D6 labels remain relevant limitations; this result does not
identify the exact hidden-label or model mechanism causing the discrepancy.
Do not treat the known-label comparison as validated prospective model selection.
The candidate is preserved as informative evidence, not recommended as an
improvement over the historical best after this result. No rollback or further
upload is made, and no molecule-level hidden labels are reconstructed.

Retained current API snapshots:

- MA: `20261002T125228445981Z-23dfbeb7`.
- D6: `20261002T125228445980Z-10056d8d`.
- A4: `20261002T125228445972Z-6fdebba0`.

All reside under `artifacts/api/classification/`. Timestamp-matched before
snapshots are `20261002T120826933682Z-e7ae47bc`,
`20261002T120831202383Z-37f7301d`, and `20261002T120835483136Z-ffb79918`,
respectively, declared in the attempt receipt. Full boards and private profile
data remain local; this report preserves the participant's shareable result.
Next classification eligibility remains October 3 at 09:08:22.746783 JST.
