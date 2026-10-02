# Classification AP-trunk probe: severe D6 Live loss, exact A4 restoration

Scored submission `20261001T141041310215Z-c7faa79f` at 2026-10-01 14:13 UTC.

| Endpoint | Metric | Previous | Current | Delta |
|---|---|---:|---:|---:|
| MA | MCC | 0.396471 | 0.317448 | -0.079023 |
| MA | Accuracy | 0.834978 | 0.833447 | -0.001531 |
| MA | Precision | 0.458982 | 0.393644 | -0.065338 |
| MA | Recall | 0.561073 | 0.460800 | -0.100273 |
| D6 | MCC | 0.308635 | 0.104986 | -0.203649 |
| D6 | Accuracy | 0.904183 | 0.871538 | -0.032645 |
| D6 | Precision | 0.341206 | 0.158841 | -0.182365 |
| D6 | Recall | 0.383490 | 0.192393 | -0.191098 |
| A4 | MCC | 0.484306 | 0.529911 | +0.045604 |
| A4 | Accuracy | 0.765773 | 0.795355 | +0.029582 |
| A4 | Precision | 0.576758 | 0.628446 | +0.051689 |
| A4 | Recall | 0.738656 | 0.729207 | -0.009449 |

D6 keeps 49 blind positives but swaps 21 each way. Its ranking AP gain in conditional labelled-only OOF did not transfer to useful Live calls. The earlier historical-rate OOF MCC point estimate was slightly negative; it did not resolve the magnitude of this prospective loss. No hidden molecule labels are inferred.

A4 is bitwise restored to the prior best 249-call vector and exactly reproduces its prior MCC/accuracy/precision/recall. This supports the intended restoration and continuity for A4, not an independent D6 generalization result.

Reject this D6 composition for current Live. Preserve it as an informative negative result and retain the original D6 vector as the performance control. Any rollback is a separate future submission requiring user approval and cooldown. No new upload, training, extraction or regression change in this audit.
