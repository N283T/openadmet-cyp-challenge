# CYP work consolidation through October 2, 2026

## Changes retained

- Structure/trunk work: four CYP reference complexes, HEM-aware multi-block
  inputs, frozen extraction provenance, pooling/readouts and incumbent-matched
  comparisons. Large external inputs remain local; manageable project-authored
  results are under `artifacts/tracked/`.
- Regression: adopted trunk replacements and scored submission records.
  [Current composition](../../configs/regression/current-boltz-adoption.json)
  improves Live 1A2 ST-RAE to .506911 and D6 to .535721, with macro .458311.
- Classification: preserved the failed AP-trunk blend, independent D6 fits,
  primary-selection/masking controls, nested threshold experiments, ten saved
  CheMeleon readouts, complete eleven-atom refits and scored submissions.
  Failures and unsent candidates remain explicit evidence.
- Infrastructure: recorded immutable hardlink input mode with hash checks,
  failed-run tombstones, associated tests and reproducible launchers/configs.

## Classification state at consolidation

| Item | State |
| --- | --- |
| Latest sent classification payload | `20261002T092348406577Z-4ea606c1`, scored |
| D6 submitted calls / Live MCC | 175 / .117740060 |
| D6 historical best MCC | .308635; no new restoration upload |
| A4 submitted calls / Live MCC | 249 / .529910537, rank 1 |
| Latest MA MCC | .323825299 |
| Actual new-model top-49 diagnostic | 22 shared with best; 27 additions and 27 removals |

The [scored result](2026-10-02-d6-chemeleon-live-result.md),
[decision audit](2026-10-02-d6-validator-decision-audit.md) and
[actual blind overlap](2026-10-02-d6-chemeleon-top49-overlap.md) supersede the
pre-submission recommendation. Known-label CV superiority did not establish
prospective Live superiority. Historical discordance, decision-policy confounding
and unobserved weak-primary outcomes remain unresolved. Neither the old model
nor 49 calls is adopted as a new modeling constraint.

## Preservation and checks

The [research log](../research/README.md) indexes hypotheses, run IDs, results
and decisions, including RL0206 through RL0211. Completed compact outputs and
own weights retain input/config/output checksums. The interrupted legacy run's
completed scientific block is included in the retained merged result; its
original local record remains intact. Private submission profiles, complete
leaderboard downloads, large run snapshots and third-party model assets stay
outside tracked files. No experimental records are deleted during consolidation.

Validation: `pixi run check` passes; `pixi run test` reports 73 passed and eight
skipped database integration tests. The repository remains private. This
consolidation changes repository organization and Git state; it makes no further
challenge submission or model-selection claim.
