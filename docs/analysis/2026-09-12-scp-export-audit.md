# Independent audit of recovered native Chemprop predictions

The first training adapter used bare target names (`CYP2D6`) when selecting
canonical labels (`CYP2D6_pIC50_direct_inhibition`). Consequently its original
OOF exports were empty and test exports used bare names. This was an export
defect: the native model's named all-head validation and test arrays were saved.
The separate recovery runs remap those immutable arrays without retraining.

Independent script `scripts/scp_repro_export_audit.py` uses literal canonical
target names and hard total counts, not a potentially empty expected-label
selection. It checks canonical compound/endpoint/fold membership, label and CI
identity, raw validation membership, exact named-head prediction mapping, all
750 official test identities per endpoint, all five fold predictions, uniqueness,
finiteness and the recomputed five-fold test mean. All inputs are declared in
the recorded runs via `scripts/scp_repro_launch_audit.py`.

| Training | Recovery | Independent audit | OOF rows | Test rows | Native array mapping |
|---|---|---|---:|---:|---|
| `fec30715`, parent single | `0f300d9d` | `b1080b38` | 1,457 | 750 | Exact |
| `03b86209`, source-scaffold single | `ea687d29` | `edea60d1` | 1,493 | 750 | Exact |
| `423d2caf`, parent auxiliary8 | `dcf54932` | `6b3569b5` | 6,384 | 3,000 | Exact |
| `e6e93fbf`, parent isoform6 | `e6fbb180` | `e5654aa2` | 1,457 | 750 | Exact |
| `288cae8b`, parent union18 | `d5c85419` | `08b26a21` | 6,384 | 3,000 | Exact |

Completed audits also reproduce the test mean within 1e-12. These checks validate
the published-array-to-canonical-export mapping; they
do not claim an independent rerun of neural inference from each checkpoint.
The source-scaffold control deliberately retains the author's grouping behavior;
the separate data/membership audit records its two cross-fold connectivity groups.

Initial audit attempts failed at the final finiteness check because NumPy
`isfinite` does not accept this Polars Series representation. The corrected script
uses Polars `is_finite`; failed runs remain recorded, including `ec6b73c1`.
This audit-code exception did not indicate a nonfinite model prediction.
