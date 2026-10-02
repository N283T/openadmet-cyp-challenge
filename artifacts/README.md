# Run artifacts and tracked deliverables

The run store is ignored by Git and immutable while retained. Completed runs,
submission records and scientifically informative negative results are durable.
Failed or interrupted runs with no usable result may be removed after their reason
is recorded; cleanup must be targeted rather than a broad deletion of ignored files.

- `runs/<track>/<id>/`: intent, source/environment/input snapshots, manifest, console,
  outputs, persistent scratch, result (including failure).
- `pruned-runs/<track>/<id>.json`: compact tombstones retained after failed run payloads
  are pruned with `pixi run cyp prune-runs --status failed --execute`.
- `boltz-trunk/post-r1/<target>/`: canonical per-record post-trunk recycle-one tensors
  and keyed Parquet manifests; legacy run output paths are non-duplicating symlinks.
- `tracked/cyp-boltz-*`: compact pooling, incumbent comparison, diversity, adoption
  and Live audit outputs. The canonical final incumbent comparison ends in
  `04ff278c`; earlier comparison snapshots are retained provenance, not separate
  current model recommendations. See [current adoption](../docs/analysis/2026-10-01-boltz-regression-adoption.md).
- `candidates/<track>/<id>/`: validated frozen prediction file and source run IDs.
- `submissions/<track>/<id>/`: separate frozen upload file, initial manifest and events.
  A prepared file is not a sent or accepted submission.
- `api/<id>/`: API contract/source snapshots, response status, revision and hashes.
- `data-imports/<id>/`: input inventories and provenance.
- `backups/<id>/`: PostgreSQL custom-format dump and checksums.
- `tracked/`: finalized project-created weights, compact derived data and runtime
  results that are safe, legally shareable and small enough for ordinary Git.

Do not place secrets, personal information, proprietary inputs, or third-party
copyrighted model/data payloads in `tracked/`. Do not use Releases or Git LFS merely
to bypass normal repository size limits. Oversized local artifacts should retain
provenance and checksums in a tracked manifest or report.

These local copies survive `/tmp` cleanup but not disk failure. Back up this tree,
`data/`, and source control to a separate disk/service before expensive campaigns.
