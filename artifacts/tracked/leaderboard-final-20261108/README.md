# Final leaderboard snapshot (2026-11-08)

The 24 public leaderboard CSVs in `csv/` are byte-for-byte copies of files fetched
from the OpenADMET Hugging Face Space. They cover both tracks, the overall board,
and every endpoint for three boards:

- `final/`: the Final Leaderboard, captured 2026-11-08 UTC after the Space moved to
  phase 4 ("challenge closed, final leaderboard released");
- `interim/`: the frozen Interim Leaderboard, captured at the same time;
- `live/`: the last Live Leaderboard, captured 2026-11-03 06:53 UTC, after the final
  submissions of both tracks. The Live tab no longer exists in phase 4.

[`manifest.json`](manifest.json) preserves each API call, capture time, original
local API path, SHA-256, byte count and row count. [`tracked-files.json`](tracked-files.json)
maps those original paths to the portable copies. Each copy was checked against
the original SHA-256 and byte count by `scripts/promote_leaderboard_snapshot.py`.
The ignored `artifacts/api/` tree retains the full API response and configuration
snapshots, including the phase-4 Space source in
`artifacts/api/20261108T123630779301Z-5464cfae/`. No leaderboard values were edited.

The result and its interpretation are in `docs/research/2026-10.md` (RL0239).
