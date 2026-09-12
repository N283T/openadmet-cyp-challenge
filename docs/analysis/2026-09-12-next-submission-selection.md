# Next submission: all-endpoint current-best request

The user authorized submission, then explicitly corrected the scope: submit the
best current configuration across all four endpoints, not the D6-only placement
control. No prediction file was uploaded during this request. The existing
D6-only and old-five-mean/all-placement files must not be silently described as
the latest global ensemble.

## Current rate limit

Fresh official source/API snapshot `artifacts/api/20260912T113927495541Z-f38c3640`
completed successfully. `config.py` defines `HOURS_BETWEEN_SUBMISSIONS = 12`.
`submission.py` independently checks the most recent submission date per track
and user before accepting another submission. The last local confirmed send is
`20260912T074820067880Z-3d2e5c50`, started2026-09-12 07:48:20UTC. Thus the
next regression eligibility is approximately2026-09-12 19:48UTC, or September13
04:48JST;04:49JST provides a margin for the receipt time. This is not a promise
of server-side eligibility if another submission is made outside this task.
No cooldown attempt or retry was made. No automatic future upload is scheduled.

## Deployment distinction

The old five-model full-training mean was sent earlier. Global126 (`84c28266`)
and174 (`7f98ba4d`) comparisons were subsequently built from saved OOF vectors;
they do not have a completed aggregate test artifact. The later native mean and
its fitted nonlinear calibration variants do have complete test predictions.
Using an old submitted vector merely because its file exists would not implement
the user's corrected scope.

A bounded deployment inventory is checking which selected global components
have complete feature coverage and saved fold models, and which need further
inference or refitting. At least the geometry mean/std, joint geometry/substrate
and focus primary components have saved fold models; absence of aggregate test
predictions does not mean their trained models are missing. Full-data refitting
must be distinguished explicitly from five-fold prediction averaging.

Selection should consider all endpoint results and apply the documented blind
placement assumptions separately from local label calibration. A best screened
OOF score is not a verified best leaderboard score. Freeze a new audited,
all-endpoint submission artifact once deployment is complete; do not relabel the
old control file. The user has already authorized sending the intended final
all-endpoint artifact, subject to actual submission availability.

Recorded inventory `20260912T114524817386Z-0456dff2` now preserves41 weighted
endpoint supports /30 selected candidate IDs and36 recursively referenced runs.
It records exact source input paths/hashes and saved asset filenames. Follow-up
inspection narrows the initial missing-baseline concern: the audit-derived
baseline has upstream saved models, so a missing model in an audit directory is
not proof refitting is required. Training-only nested auxiliary tables still
need test routing from the underlying all-molecule predictions. The aggregate
all-endpoint test artifact remains incomplete; no candidate/submission has been
registered for this new scope.
