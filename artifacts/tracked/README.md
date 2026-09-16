# Tracked artifacts

This directory holds finalized project-created artifacts that are appropriate for
ordinary Git: compact derived datasets, model weights and runtime results.

Before adding a file, verify that it:

- contains no secret or personal information;
- is project-authored or legally redistributable;
- does not embed a third-party copyrighted model or dataset;
- fits normal Git hosting limits without Git LFS or a Release workaround; and
- has clear provenance in a report, config, or manifest.

Transient run snapshots remain under the ignored `artifacts/runs/` tree. Promote only
the finalized files that are useful to collaborators.
