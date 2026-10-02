# Exact new-model top-49 overlap with the historical best

The user asks how much the submitted model would differ from the best-scored
49-call vector if its decision were limited to its top 49 molecules. This is
an identity comparison on the actual 750 blinded molecules, not an OOF estimate.

Recorded run `20261002T132058268073Z-dfc673b6` ranks the saved ten-model mean
margin, verifies that ranking by mean raw score produces the identical top set,
and exactly reproduces every actually submitted 175-call decision. Molecule
IDs/names match both immutable payloads. The historical best reference is
prepared `20260924T050705396203Z-10e1d3ad` (49 D6 positives, scored MCC .308635).
The current submitted reference is `20261002T092348406577Z-4ea606c1`.

| New-model decision | Positives | Shared with best | Added | Removed | Changed decisions |
| --- | ---: | ---: | ---: | ---: | ---: |
| Counterfactual top 49 | 49 | 22 | 27 | 27 | 54 |
| Actual submitted policy | 175 | 43 | 132 | 6 | 138 |

The top-49 set retains only 44.9% of the historical best's positives. Of the
27 displaced historical positives, 21 rank between 50 and 175 in the new
model, while six fall beyond 175. A score/margin gap exists between ranks 49
and 50, so tie handling does not affect the result.

The 175-call candidate expands the set substantially, while its narrow top set
also changes substantially. These are separate observed properties. They do
not identify which additions/removals are errors or what an unsubmitted top-49
payload's Live MCC would be. The top-49 vector is a diagnostic, not an approved
submission or new fixed positive budget. No model is fitted or upload made.

Retained keyed comparisons, all 54 changed identities, and all historical-best
new ranks are under
`artifacts/tracked/cyp-d6-zero-start-20261002T132058268073Z-dfc673b6/`.
