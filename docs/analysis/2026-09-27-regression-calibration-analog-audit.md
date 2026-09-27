# Why the affine calibration OOF gains did not transfer

**2026-09-27 (JST).** The September 27 Live regression submission showed that
the 1A2 and 3A4 affine changes worsened ST-RAE despite improvements on the
full parent OOF. This follow-up stays within calibration: it checks whether
the submitted transformations behave differently on an existing subset of
training analogs chosen around strong 1A2/2C9/3A4 hits. This is closer to the
official test acquisition than the full labelled training population, but it
does not establish the blind scored population or its optimal placement.

Run `regression/20260927T054145210148Z-ba3de34b` applies each fixed affine
map around the complete OOF prediction mean, exactly as in the earlier
screens, and scores the same models on all labelled parent OOF rows and on
the analog subset. It uses 198 of 1,372 1A2 rows and 238 of 2,300 3A4 rows.
The subset mean pIC50 is higher than the full set for both targets:
`5.129` versus `4.908` for 1A2 and `4.403` versus `4.066` for 3A4.

| Target and fixed map | Full OOF ST-RAE change | Analog OOF change | Live change |
| --- | ---: | ---: | ---: |
| 1A2, scale `.925` | −.001121 | **+.003707** | **+.003815** |
| 3A4, shift `+.05`, scale `.9625` | −.006659 | **+.003423** | **+.000435** |
| 3A4, shift `+.05` only | −.005485 | −.000162 | Unmeasured alone |
| 3A4, scale `.9625` only | +.001045 | **+.006938** | Unmeasured alone |

The analog subset correctly gives the *direction* of both submitted losses,
unlike the full OOF. In its five folds, 3A4 narrowing alone worsens ST-RAE
in **5/5**; a 1,000-replicate molecular-series bootstrap for this fixed map
has a conditional 95% delta interval `[+.00088,+.01282]`. The shift-only 3A4
effect is effectively zero on the analog subset, with interval
`[−.01097,+.01153]`. Thus narrowing is a plausible source of the 3A4 loss,
but the Live submission changed center and spread together, so its separate
effects are not identified by that result. For 1A2, narrowing worsens three
of five analog folds and its interval crosses zero (`[−.00660,+.01257]`).
A `.075` widening of 1A2 is also slightly worse on the analog subset
(`+.000847`), providing no supported correction direction from this audit.

The 75 hit anchors were selected using 1A2/2C9/3A4 potency; analog labels
are sparse, and the subset may not match the actual blind scored analogs.
The bootstrap conditions on this subset, saved model predictions, and maps;
it does not account for proxy construction or the earlier calibration search.
These results are a **diagnostic guardrail**, not evidence that another affine
change will improve Live ST-RAE. Retain the September 26 regression payload
as performance control. Calibration remains the focus, but a further upload
should isolate a single center or spread change and have a stronger
test-relevant rationale than the full-OOF grid optimum alone.
