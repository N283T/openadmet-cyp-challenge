# Run artifacts and tracked deliverables

The run store is ignored by Git and immutable while retained. Completed runs,
submission records and scientifically informative negative results are durable.
Failed or interrupted runs with no usable result may be removed after their reason
is recorded; cleanup must be targeted rather than a broad deletion of ignored files.

- `runs/<track>/<id>/`: intent, source/environment/input snapshots, manifest, console,
  outputs, persistent scratch, result (including failure).
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
