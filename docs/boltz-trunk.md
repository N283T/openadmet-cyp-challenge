# Boltz trunk input configuration

The primary CYP trunk experiment is a template-free Boltz-2 encoder run, not
an affinity prediction. Each record contains the same fixed receptor context
and one variable challenge compound:

- protein chain `A`: full canonical CYP sequence with the pinned LiteFold A3M;
- ligand chain `H`: fixed CCD component `HEM`;
- ligand chain `B`: one standardized challenge SMILES; and
- one soft pocket constraint from `B` to the target's extended contact-residue
  set (`max_distance=8.0`, `force=false`).

There is no `properties` section and no template section in the generated
Boltz YAML. Therefore neither the affinity head nor a structural template is
requested. The structure checkpoint `boltz2_conf.ckpt` supplies the frozen
trunk. The optional template comparison is deferred to a separate configuration
so it cannot silently alter the primary feature definition.

Runtime v1 uses three recycling steps, CUDA kernels, BF16 mixed precision,
FlashAttention disabled, a maximum of 8,192 parsed MSA rows, and deterministic
1,024-row MSA subsampling per pass with seed `20260928`. RDKit ETKDGv3/UFF
ligand conformers use the same explicit seed. The fixed HEM is a
normal non-polymer entity, not the pocket binder; the variable compound `B` is
the only pocket binder.
Only post-trunk `s` and `z` tensors are retained, in float16. Pre-trunk tensors
are not written. The extractor exposes `pre`, `post`, and `both` stage choices
for later ablations without changing the primary v1 definition.

Completed recycle-one outputs are consolidated under
`artifacts/boltz-trunk/post-r1/<target>/records/`. The root and target-specific
`manifest.parquet` files map record filenames to compound IDs, roles, standardized
SMILES, source runs or archives, byte sizes and SHA-256 hashes. Legacy run output
paths are relative symlinks to this canonical store and do not duplicate tensors.

The immutable settings are in
[`configs/boltz-trunk/prepare-v1.json`](../configs/boltz-trunk/prepare-v1.json).
`scripts/launch_boltz_trunk_input_prep.py` records expansion of the fixed
context and the target-specific compound set into one standard Boltz YAML per
compound. Each CYP uses only compounds with at least one target-specific train
measurement plus all 750 test compounds: 2,163 / 2,035 / 2,247 / 4,334 records
for CYP1A2 / CYP2C9 / CYP2D6 / CYP3A4 respectively.
The current full input runs are `shared/20260929T034044572251Z-b2e9435f`,
`shared/20260929T034046886086Z-b2f6b507`,
`shared/20260929T034049410496Z-672bc9b5`, and
`shared/20260929T034052097919Z-b106051b` in CYP1A2 / CYP2C9 / CYP2D6 /
CYP3A4 order. Their manifests contain no duplicate compound/role pairs and no
missing YAML records.

CYP3A4 GPU smoke run `shared/20260929T033909092009Z-3b347c28` completed one
train and one test record through the structure checkpoint. It used the
cuequivariance CUDA triangle kernels, retained only post-trunk `s` and `z`, and
wrote both tensors as float16. The two compressed files were approximately
75 MiB each. A recycle=3 CYP2C9 extraction was stopped after 870 records when
the primary configuration was changed to recycle=1; its run directory was
deleted by explicit decision. Recycle=1 smoke run
`shared/20260929T062141363586Z-52850879` completed before restarting the full
extraction. All four target extractions are now complete in the canonical post-r1 store.

## Downstream pooling comparison (October 1, 2026 JST)

The canonical manifest now contains all 10,779 records. The
[fixed experimental protocol](analysis/2026-10-01-boltz-trunk-pooling-protocol.md)
compares 12 PXR-style pooled feature views with separate/mixed HEM and Fe
blocks, trunk-only and molecular-augmented features, LGBM and TabPFN, and all
five UMAP folds. Pooling and fitting are recorded runs. Large retained input
tensors use explicit hash-checked hardlink snapshots; source tensors remain
immutable. TDI results are diagnostic because the complete operating pool lacks
trunk features. No submission is part of this comparison.

The UMAP regression screen is exploratory. The subsequent
[incumbent-matched comparison](analysis/2026-10-01-boltz-incumbent-redo.md)
uses exact submitted parent-series cohorts, labels and folds. The final comparison
run is `regression/20261001T111723965954Z-04ff278c`; earlier report snapshots
remain retained for provenance. [Trunk diversity and multiple-member placement](../artifacts/tracked/cyp-boltz-multi-report-20261001T113928351330Z-e5c62212/report.md)
follow that same protocol.

The user subsequently authorized regression adoption and submission. The
[adoption and Live result](analysis/2026-10-01-boltz-regression-adoption.md)
documents the submitted 1A2/D6 composition, exact C9/A4 rollback and prospective
score improvements. This is the current regression checkpoint. Classification
deployment is the next separate task; incomplete-pool TDI diagnostics do not
justify a positive-call rate or submission composition.
