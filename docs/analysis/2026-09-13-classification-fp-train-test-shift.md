# Train–test Morgan fingerprint coverage

## Conclusion

The blind test molecules are slightly less fingerprint-dense than the labeled
sets, while the hashed 2,048-bit position vocabulary is almost completely covered
by training. That apparent coverage is misleading: unhashed Morgan identifiers
show substantial test novelty, especially against the small CYP2D6 set. D6 also
has much weaker nearest-neighbor support than A4. The issue is therefore not a
fingerprint-generation mismatch or missing region of the 2,048-bit array; it is
frequency shift plus unresolved structural environments hidden by hashing.

## Design

Diagnostic run `20260913T131749408718Z-0167f830` compares the 750 blind test
molecules with the 1,497 D6 labels, 3,584 A4 labels and all 6,146 unblinded
molecules. No test labels are used. Measurements include:

- Active hashed bit positions and per-molecule density.
- Radius-0, radius-1 and radius-2 occurrence density.
- Per-position prevalence correlation, Jensen-Shannon divergence and largest
  shifts.
- Unhashed sparse Morgan environment-ID novelty, both unique and
  occurrence-weighted.
- Maximum test-to-train Tanimoto similarity.

## Density and hashed positions

| Train reference | Mean active bits | Test mean | Density change | Train/test active positions | Test-only positions |
|---|---:|---:|---:|---:|---:|
| All unblinded, n=6,146 | 49.09 | 46.23 | -5.8% | 2,048 / 1,774 | 0 |
| D6 labels, n=1,497 | 48.04 | 46.23 | -3.8% | 2,039 / 1,774 | 6 |
| A4 labels, n=3,584 | 49.80 | 46.23 | -7.2% | 2,046 / 1,774 | 1 |

The six D6-test-only positions account for only 0.063% of test on-bit
occurrences; the single A4-test-only position accounts for 0.003%. Broad 128-bit
index-bin occurrence mass differs by at most about 1.3 percentage points against
all unblinded training. This argues against a bit-ordering, fingerprint version or
pipeline-alignment failure.

Test density is lower mostly beyond the center atom. Against D6, the mean test
minus train changes are -0.07 radius-0 bits, -0.76 radius-1 bits and -1.03
radius-2 bits per molecule. Against A4 they are -0.36, -1.64 and -1.63. Test
molecules therefore contain fewer recorded local and extended environments,
rather than a markedly different number of atom-identity bits alone.

## Frequency shift

Per-bit train/test prevalence rank correlation is 0.686 for both endpoint sets.
There are 209 D6 bit positions with absolute prevalence shift above two percentage
points and 46 above five points. A4 has 225 and 59 respectively. Occurrence-based
Jensen-Shannon divergence is 0.072 bits for D6 and 0.077 for A4.

Several large shifts affect bits already implicated by the attribution analysis.
The two FDR-significant negative D6 bits become more common in test: bit 1654
rises from 7.0% to 12.1%, and bit 1872 from 3.7% to 8.3%. For A4, positively
associated bit 1114 rises from 20.7% to 27.2%, while strongly negative bit 294
falls from 21.2% to 6.9%. These directions could alter probability calibration
and ranking, although the fixed-count submission policy prevents them from
directly changing the number of positive calls.

## Novel unhashed environments and neighborhood support

| Train reference | Test sparse IDs unseen in train | Test occurrences from unseen IDs | Median nearest Tanimoto | Test below 0.4 |
|---|---:|---:|---:|---:|
| All unblinded | 16.8% | 1.53% | 0.587 | 0.13% |
| D6 labels | 34.3% | 4.58% | 0.471 | 42.0% |
| A4 labels | 26.6% | 3.03% | 0.540 | 18.9% |

Unique-ID novelty is high because rare environments form a long tail; the
occurrence-weighted numbers are the more relevant magnitude. Even by that measure,
D6 has roughly three times the novelty of the full unblinded pool. Its Tanimoto
coverage is also materially worse: 12.5% of test molecules are below 0.3 maximum
similarity, compared with 2.5% for A4.

## Interpretation

The fixed 2,048-bit representation makes train and test look position-compatible,
but collisions fold genuinely unseen environments into familiar locations. The
D6 weakness is consistent with small labeled-set coverage rather than excessive
raw bit sparsity. Useful follow-ups are larger or unhashed fingerprints, explicit
rarity/support features, and blending a local similarity model with a global
model. Density alone should not be rescaled, and the unsupervised test comparison
does not establish which shifted environments are TDI-positive.
