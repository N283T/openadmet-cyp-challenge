# Experiments derived from the CYP research reports

The user authorized testing the four remaining directions, with regression and
CYP1A2/CYP2D6 as priorities. Existing CYP2C9/CYP3A4 and auxiliary labels remain
usable when helpful. Every comparison retains the existing parent split and
records all inputs, successful/failed attempts, fitted models and keyed OOFs.

## Prespecified program

1. **TDI continuous supervision.** Add measured TDI-condition potency, then its
   observed shift from direct potency, to the existing external-pretrained
   multitask encoder. Match the original interval-loss training recipe and use
   identical new-primary LGBM readouts. See [the TDI experiment](2026-09-12-report-tdi.md).
2. **Constrained response model.** Compare independent shared-MLP heads with
   condition-matched Hill-linked response heads. The primary screen used the
   active TDI arm, so the relation uses TDI potency and TDI Emax, while direct
   pIC50 remains the main supervised endpoint. Use frozen external600+47 inputs,
   the same grouped inner epoch selection and outer refit, and fixed slope versus
   a bounded per-CYP learned slope. Per-CYP affine scale accounts approximately
   for the screen's processed log2fc scale. This is an observation surrogate,
   not a claim that public log2fc equals a raw inhibition fraction.
3. **Ligand geometry and charge.** Compare explicit topology/charge proxies and
   ETKDG geometric descriptors beyond Jazzy, using matched LGBM and TabPFN.
   No pKa or pH-microstate prediction is claimed. See
   [the geometry experiment](2026-09-12-report-geometry.md).
4. **Substrate supervision.** Use real external substrate labels, exclude all
   challenge connectivity, and add frozen substrate probabilities. This is the
   practical substrate branch of the substrate/SoM proposal; atom-level SoM is
   not tested in this round. See [the substrate experiment](2026-09-12-report-metabolism.md).

After the individual comparisons, rerun the same global Caruana selection with
new compatible OOFs. Compare with both the previous manual candidate and the
previous global result. Weight cross-fitting remains a weight-stability check,
not newly nested base-model training. No submission is part of this round.

## Observation-data preparation

Successful preparation: `20260912T093252891365Z-271e3020`. Direct targets retain
original measurements; TDI potency and native signed TDI Emax are kept separately.
Actual molar screening concentration is preserved. Finite duplicate observations
are collapsed by compound/enzyme using a median; missing target values stay missing.

Two initial attempts (`a7129af9`, `3c238d50`) rejected an unmatched platinum compound
in the older local Emax file. The successful run uses the already verified refreshed
Emax snapshot from `fd93c411`; no molecular identity was guessed or reassigned.
Both failed records and the exact unmatched source row remain available.

For paired observed labels, fixed-slope Hill response correlates with measured
screen log2fc at Spearman 0.908/0.935/0.842/0.969 for 1A2/2C9/2D6/3A4. These are
observed-label consistency checks, not model validation scores. They motivate
testing the relation in a controlled learner.

## Completed results

All four branches above completed, plus simultaneous geometry/substrate/Jazzy
readouts. Atom-level SoM remains untested.

| Matched comparison (ST-RAE, lower is better) | CYP1A2 | CYP2D6 |
|---|---:|---:|
| TDI control LGBM | 0.692959 | 0.837797 |
| TDI arm + measured shift LGBM | 0.690031 | 0.835201 |
| Geometry control TabPFN | 0.705017 | 0.874374 |
| Geometry mean/std TabPFN | 0.692620 | 0.877186 |
| Substrate control LGBM | 0.691998 | 0.842287 |
| Six substrate scores LGBM | 0.689724 | 0.842759 |
| Independent observation MLP | 0.734246 | 0.886731 |
| Fixed-Hill observation MLP | 0.743712 | 0.904206 |
| Learned-Hill observation MLP | 0.737859 | 0.894060 |

Hill constraints were negative in this implementation; its independent control
is also weak. Observation training `ec2b7210` and independent audit `8df31bae`
verify six OOFs, all 15 checkpoint scalers, memberships and CPU reconstruction.
Selected and refitted epochs were 1–5, verified against history. See
[implementation](2026-09-12-observation-model-implementation.md) and
[audit, including failed attempts](2026-09-12-report-observation-audit.md).

Joint preparation `b780d545`, readout `4ec73c9f`, audit `9ec87ba4` combine
geometry/substrate, then Jazzy. Geometry plus substrate improves standalone
2D6 from 0.842287 to 0.837145, but fixed additions do not improve manual macro.

Expanded global selection `84c28266`, audited by `846181c4`, adds 19 candidate
identities (107 to 126). Weight-crossfit macro moves 0.583257 to 0.582816;
1A2 moves 0.674078 to 0.670794 and 2D6 0.825782 to 0.825436. Other endpoints
slightly worsen. Full-OOF fitted weights score 0.572518, separately reported
because the same labels select weights. See
[aggregate comparison](2026-09-12-report-global-ensemble.md).

Retain all assets. No substantial focus-target breakthrough or new submission.
