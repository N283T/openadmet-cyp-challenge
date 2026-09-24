# Classification candidate preflight

This is a retrospective diagnostic for a fixed-call candidate. The current
[candidate evaluation protocol](classification-evaluation.md) selects both the
model and positive-call rate; this preflight alone cannot make that choice. It
compares a fixed-call candidate with the current control using only training
labels, honest OOF predictions and blind-test structures. It does not use hidden
test labels or leaderboard results.

## Required inputs

- `decision-test.parquet`: long keyed control/candidate scores and calls with
  `compound_id`, `model`, `score`, `call`;
- `decision-oof.parquet`: the same fields plus `y_true`;
- `member-test.parquet` and `member-oof.parquet`: long keyed ensemble-member
  predictions and calls with `compound_id`, `model`, `prediction`, `call`;
- canonical `compounds.parquet` and `measurements.parquet`; and
- a config naming the endpoint, control, candidate, member models, call count
  and thresholds.

Launch through the recorded-run wrapper:

```sh
pixi run python scripts/classification_candidate_preflight_launch.py \
  --config configs/analysis/MY-PREFLIGHT.json \
  --decision-test PATH/decision-test.parquet \
  --decision-oof PATH/decision-oof.parquet \
  --member-test PATH/member-test.parquet \
  --member-oof PATH/member-oof.parquet \
  --compounds PATH/compounds.parquet \
  --measurements PATH/measurements.parquet
```

Use explicit input paths from completed runs. The wrapper snapshots every input
through `cyp run` and writes `summary.json`, `checks.csv`, pairwise member
dependence, OOF swap yields and blind swapped-neighborhood diagnostics.

## Default gates

The initial conservative defaults are:

| Gate | Default |
|---|---:|
| Fraction of current positive calls replaced | <= 25% |
| Mean blind-test member Spearman | <= .95 |
| Effective independent vote count | >= 2.0 |
| OOF swap positive-rate advantage | >= .02 |
| OOF swap Fisher two-sided p | <= .10 |
| Test-minus-OOF member Spearman increase | <= .05 |
| Removed-control label-neighborhood advantage | veto at p <= .05 |

Set `"enforce": true` only when explicitly using these conservative checks as a
required stop for a fixed-call comparison. Any failed gate exits nonzero after
writing diagnostic outputs. Use `"enforce": false` for descriptive candidate
evaluation or retrospective replay.

A failed gate is not proof that the candidate lacks value. These historical
thresholds were not calibrated to maximize leaderboard gain or to select a call
rate. Record any decision to use or override them in the research log. Passing
preflight likewise does not predict leaderboard improvement.

## Retrospective validation

Run `20260918T104912578742Z-11a49571` replays the preflight against the failed
D6 support6 submission. It fails all six checks: swap fraction `.551`, member
Spearman `.970`, effective votes `1.13`, OOF swap Fisher `p=.865`, correlation
inflation `.114`, and label-neighborhood veto `p=.0005`. The preflight would
therefore have blocked candidate promotion before upload.
