# September 26 leaderboard snapshot

The 16 public Live and Interim leaderboard CSVs in `csv/` are byte-for-byte
copies of the files fetched from the OpenADMET Hugging Face Space on
2026-09-26 UTC. They cover both tracks, the overall board, and every endpoint.

[`manifest.json`](manifest.json) preserves each API call, capture time, original
local API path, SHA-256, and row count. [`tracked-files.json`](tracked-files.json)
maps those original paths to the portable copies. Each copy was checked against
the original SHA-256 and byte count before being added to Git. The ignored
`artifacts/api/` tree retains the full local API response and configuration
snapshots. No leaderboard values were edited in the tracked copies.

Use the paths in `tracked-files.json` when rerunning the analysis from a fresh
checkout. The analysis reports are in `docs/analysis/2026-09-26-*` and the
decisions are in `docs/research/2026-09.md` (RL0168–RL0171).
