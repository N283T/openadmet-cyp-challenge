# Independent CYP2D6 measurement coverage near the known hit anchors

The recent D6 weak-activity interventions did not improve ST-RAE on both
splits. This audit asks a narrower question: whether the existing external
corpus contains enough independent, nearby D6 measurements to validate a
different weak-activity strategy before further model work.

## Method

Recorded run `regression/20260926T030151905020Z-d1a709d3` uses the 70 unique
challenge hit anchors in the parent assignments. For each external compound
with a D6 observation, it computes the maximum Morgan radius-2, 2048-bit
Tanimoto similarity to those anchors. Challenge-overlap compounds are
excluded. Coverage is counted separately by source and measurement type;
exact values and upper-censored values are kept distinct. Weak means an
exact p-scale value below 4 or a censored upper bound at most 4. Similarity
to an anchor is a chemical proximity screen, not proof of identical assay
behavior or of membership in the hidden test distribution.

The four sources span 23,664 unique measured compounds and 30,464 D6
observations before the similarity cutoff. ChEMBL contains multiple assays;
each of the other sources is one assay in this corpus. IC50, Ki, pIC50 and
AC50 are reported separately and must not be pooled as equivalent labels.

## Coverage

| Source and endpoint | Similarity ≥0.3: compounds / anchors | Similarity ≥0.4: compounds / assays / anchors | ≥0.4 exact weak / censored upper ≤4 | ≥0.5 compounds |
| --- | ---: | ---: | ---: | ---: |
| ChEMBL36 IC50 | 258 / 53 | 15 / 10 / 6 | 0 / 1 | 1 |
| ChEMBL36 Ki | 4 / 2 | 0 / 0 / 0 | 0 / 0 | 0 |
| Novartis pIC50 | 639 / 57 | 43 / 1 / 16 | 0 / 0 | 6 |
| PubChem AID1851 AC50 | 346 / 51 | 25 / 1 / 15 | 0 / 0 | 3 |
| Tox21 P450-Glo AC50 | 46 / 20 | 6 / 1 / 3 | 0 / 0 | 0 |

At the looser 0.3 cutoff, ChEMBL IC50 has 15 compounds with an upper bound
at most 4 but no exact weak values; PubChem has one exact weak compound.
Novartis has no exact pIC50 below 4 even over its full 15,849-compound D6
set (minimum exact value 4.033), so its lack of weak observations near the
anchors partly reflects the source's value range. Tox21 also supplies no
exact value below 4. Counts across sources may overlap and should not be
summed into a unique-compound total.

## Decision

This corpus does not provide a dense, independent weak-D6 validation set
near the 70 anchors. It cannot resolve whether a new weak-activity remedy
will generalize to challenge-like compounds. This does not rule out external
data as pretraining input or broader-domain validation data. The next useful
step for this specific failure mode is new, relevant D6 measurements or a
substantially different learning signal with an independent evaluation set.
Stop incremental ensemble-weight and scalar-correction screens for now;
retain the submitted D6 vector. No candidate or submission follows from
this audit.

Compact source-specific counts are in
[`artifacts/tracked/regression-d6-external-anchor-coverage-20260926/`](../../artifacts/tracked/regression-d6-external-anchor-coverage-20260926/README.md).
The row-level annotated external observations stay in the recorded run.
Earlier runs `regression/20260926T025908414965Z-e311afc0` (exact values only)
and `regression/20260926T025957306693Z-10ad7122` (without Tox21) document
the sequential audit; the run above is the complete four-source result.
