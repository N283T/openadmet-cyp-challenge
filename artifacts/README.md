# Durable local artifacts

Directories are append-only and ignored by Git. No automatic cleanup.

- `runs/<track>/<id>/`: intent, source/environment/input snapshots, manifest, console,
  outputs, persistent scratch, result (including failure).
- `candidates/<track>/<id>/`: validated frozen prediction file and source run IDs.
- `submissions/<track>/<id>/`: separate frozen upload file, initial manifest and events.
  A prepared file is not a sent or accepted submission.
- `api/<id>/`: API contract/source snapshots, response status, revision and hashes.
- `data-imports/<id>/`: input inventories and provenance.
- `backups/<id>/`: PostgreSQL custom-format dump and checksums.

These local copies survive `/tmp` cleanup but not disk failure. Back up this tree,
`data/`, and source control to a separate disk/service before expensive campaigns.
