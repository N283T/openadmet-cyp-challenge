# Frozen CheMeleon fingerprints

Representation: pretrained CheMeleon `BondMessagePassing` followed directly by
`MeanAggregation`, producing 2,048 float32 values per standardized compound.
No CYP labels, fine-tuning, predictor projection, scaling or dimensionality reduction
are used. The same frozen features cover both tracks, auxiliary compounds and test.
Routine downstream evaluation remains the fixed UMAP split.

The previous PXR initial 300-dimensional representation included an untrained
`RegressionFFN` projection through `MPNN.encoding()`. The later PXR raw-2048 audit
corrected this. This project follows the direct pretrained representation and
checks equivalence to Chemprop's `MPNN.fingerprint()` path, used by the
[official extractor](https://github.com/JacksonBurns/chemeleon/blob/b5e5bd88a7070c3f278e7e65e85b792617bf125c/chemeleon_fingerprint.py).

## Reproducibility

- Config: `configs/features/chemeleon.json`.
- Dependencies: PyTorch and Chemprop installed and locked through pixi.
- Weights: Zenodo record 15460715, `chemeleon_mp.pt`; existing PXR file is declared
  as a run input and SHA-256 pinned. Official checksum re-verification was unavailable
  during this extraction because Zenodo returned HTTP 504; see the recorded outcome below.
- Source reference: commit `b5e5bd88a7070c3f278e7e65e85b792617bf125c`, with original
  extractor, README and license preserved alongside the verified weight asset.
- Asset preparation: `scripts/prepare_chemeleon.py` through `pixi run cyp run`.
  It normally checks the record checksum, with a direct-download comparison fallback.
  `configs/features/chemeleon-local-asset.json` explicitly records the local-only
  provenance fallback used after both remote methods failed.
- Extraction: `python -m cyp_challenge.chemeleon` through the same wrapper; declare
  `compounds.parquet`, SHA-256-pinned checkpoint and provenance as input files.
- Seed 42, inference mode, frozen message-passing weights, deterministic PyTorch
  algorithms, TF32 disabled; runtime records package versions, device and CUDA.
- Output: keyed Parquet shards plus `embeddings.parquet`; manifest and summary
  record hashes and all feature column names. No row-index joins.

Validation covers exact compound-ID coverage, uniqueness, finite vectors and expected
dimension; eight compounds are recomputed with a different batch size, as individual
molecules, and through the official fingerprint path. Numeric tolerances and maximum
differences are recorded. This checks the extraction path, not model accuracy.

## Database

Migration `006_molecular_embeddings.sql` adds `embedding_sets` and
`compound_embeddings`, keyed by feature run and compound ID. The 2,048-dimensional
vectors are PostgreSQL `real[]`; dataset, checkpoint, dimension and validation
provenance are retained in the set metadata. Writes are transactional and append-only.

Use `pixi run python scripts/register_embeddings.py --run FEATURE_RUN_ID` after
successful extraction. The loader verifies completed-run hashes, dataset/checkpoint
input lineage and vector coverage before insertion. Loading the same run again is
idempotent. Weights and vectors remain outside Git; code, configs and this report
are tracked. Research outcomes and failures are recorded in the research log.

## Recorded extraction — 2026-09-10 UTC

- Asset preparation: `20260910T012612848389Z-bfeed434`.
- Feature run: `20260910T012630927160Z-e9576045`.
- Dataset: `20260909T132734119233Z-a7063ccd`.
- Coverage: all 6,896 standardized compounds (6,146 non-test + 750 test), 2,048 dimensions.
- Checkpoint SHA-256: `c376624d3407204e780a0ed13a9ac097cc9bb1c13ef89cdbc633c1715c183651`.
- Features SHA-256: `8f680376f5a530e8aa6dc04758946cebba1f5574a3cf9c305c35a703f99aad35`.
- PyTorch 2.14.0, Chemprop 2.2.4, CUDA 13.0, NVIDIA RTX 5080.
- All feature values finite; complete unique IDs. Rebatching / official-path maximum
  absolute difference: `1.7881393432617188e-07`; singleton difference:
  `1.1920928955078125e-07`. Both below configured tolerances.

The reusable table is
`artifacts/runs/shared/20260910T012630927160Z-e9576045/outputs/embeddings.parquet`:
`compound_id`, then `chemeleon_0000` through `chemeleon_2047`. Join by compound ID;
the dataset's aliases/test mapping supply original molecule names and submission order.
SQL: filter `compound_embeddings.feature_set_id` by the feature run above.

**Provenance limit:** Zenodo API requests and the direct file download returned HTTP
504 in three recorded preparation attempts. This extraction uses the existing PXR
checkpoint with a reproducible SHA-256, not a newly verified official download.
The successful asset manifest explicitly sets `official_checksum_verified: false`;
DB set metadata retains this provenance. The pinned official GitHub reference source
was retrieved successfully. Recheck the official checksum when Zenodo is available,
as a new recorded operation; never replace this run's weights or feature files.

See [RL-0004](research/2026-09.md#rl-0004--frozen-chemeleon-raw-2048-fingerprints)
for the failed attempts and decision to retain the local-weight extraction.
