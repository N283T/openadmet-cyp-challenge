# Molecular preprocessing and descriptors

## Scope

All five existing official CSVs (test, inhibition, TDI, Emax, single-concentration)
are inputs to a recorded shared run. This stage creates no split, fitted scaler,
imputer, feature selection, model, or predictions.

ChEMBL structure pipeline 1.2.4 standardizes each distinct raw SMILES and extracts
its parent with the library's salt/solvent lists and default neutralization rules.
Canonical isomeric SMILES preserve stereochemistry; no extra tautomer canonicalizer
is applied. Parent canonical SMILES are keyed by SHA-256. ChEMBL exclusions and
parse/standardization failures remain in raw records and a quarantine table rather
than falling back to untreated molecules and calling them standardized.

## Duplicate policy

- One descriptor calculation per unique standardized parent across all sources.
- Keep original `Molecule_Name` and raw SMILES as aliases, including all test IDs.
- Remove only exact source-record duplicates; retain every original row's mapping.
- DRC observations match on standardized compound, endpoint, value, and available
  confidence bounds/std. Exact repeats across files share a measurement ID and
  retain all source-record links. Conflicting measurements are kept separately.
- Aux rows retain enzyme, batch, plate and concentration in their original payload.
  Repeated assays are not averaged. Train/test structural overlaps are reported,
  not automatically used for fitting or removed from the submission mapping.

The initial preprocessing run is `20260909T132734119233Z-a7063ccd`:
6,897 distinct input SMILES; 6,896 usable parents; 14 changed; one excluded platinum
compound (`OCNT-2329004`, present in Emax only). No merged parent groups or test/train/
aux structural overlaps were found. There are 47,991 non-null DRC measurement
occurrences across files, reduced to 36,385 unique measurements with provenance.
No conflicting measurement groups were found. Original files remain unchanged.

## Descriptors

Compute from parent SMILES using Python libraries, not PostgreSQL descriptor
functions. The installed RDKit 2026.03.6 registers 217 descriptors in
`Descriptors._descList`; Mordred-community 2.0.7 provides 1,613 descriptors with
`ignore_3D=True`. The exact names and package versions are stored per feature run.
The RDKit cartridge stores searchable molecule objects only.

Every declared feature key is retained. Undefined values, non-finite results and
exceptions become null with a per-compound/per-descriptor reason. There is no
imputation, constant-column removal, scaling or descriptor selection at this stage.
Per-batch Parquet files persist throughout computation. Completed feature sets are
loaded atomically; incomplete runs cannot be registered as complete DB feature sets.
A failed compute run's batches remain on disk; automatic batch resume is not yet
implemented. Workers use spawn to avoid inheriting Polars thread pools through fork.

## Reproduction

The configuration is `configs/preprocessing/default.json`. Run preprocessing with
all five CSVs as separate `--input` arguments:

```sh
pixi run cyp run --track shared --purpose "Prepare official molecular data" \
  --config configs/preprocessing/default.json \
  --input data/raw/cyp-challenge-TEST-BLINDED.csv \
  --input data/raw/cyp-challenge-TRAIN_Emax.csv \
  --input data/raw/cyp-challenge-TRAIN_TDI.csv \
  --input data/raw/cyp-challenge-TRAIN_inhibition.csv \
  --input data/raw/cyp-challenge-single-concentration-TRAIN.csv \
  -- python -m cyp_challenge.preprocessing

pixi run cyp run --track shared --purpose "Compute all RDKit and Mordred 2D descriptors" \
  --config configs/preprocessing/default.json \
  --input artifacts/runs/shared/DATASET_RUN/outputs/compounds.parquet \
  -- python -m cyp_challenge.descriptors

pixi run python scripts/load_molecular_data.py \
  --dataset-run DATASET_RUN --feature-run FEATURE_RUN
```

The loader verifies the completed runs' artifact hashes and descriptor-input
lineage. It appends immutable dataset/feature versions, never truncates tables,
and is idempotent for an already-loaded run. Loader attempts are logged under
`artifacts/data-loads/` even if the transaction fails.

## Database access

- `molecular_datasets`: preprocessing versions, keyed by run ID.
- `compounds` / `dataset_compounds`: unique parents and dataset membership.
- `source_records` / `record_origins`: original payloads and row-level provenance.
- `measurements` / `measurement_sources`: deduplicated DRC values and source links.
- `single_concentration`: typed convenience view over original aux records.
- `feature_sets`: source dataset ID, descriptor catalog, versions and quality summary.
- `compound_descriptors_full` / `compound_mordred`: JSONB keyed by descriptor name,
  one row per `(feature_set_id, compound_id)`.
- `descriptor_errors`: explicit null reasons for each feature row.

Always select dataset and feature-set IDs explicitly. Repeated preprocessing runs
must not accidentally multiply training rows when joining versioned tables.

## Initial feature set and validation

Feature run: `20260909T132918829974Z-bc81af58`.
Both feature families cover all 6,896 standardized parents, including all 750 test
compounds. RDKit has zero missing values. Mordred has 1,020,722 null values;
84 columns are entirely null (atom-type extrema for atom types absent from this
collection). These columns are retained, with their missing-value reasons.
All-null/constant-column lists and the DB/Parquet reconciliation are saved in
[preprocessing-quality.json](preprocessing-quality.json).

The audit checks passed for lineage, molecule counts, descriptor columns, null
counts, all original source rows, test identities, and descriptor-error coverage.
22 local tests passed, including salt-parent handling, stereochemistry, invalid
and excluded molecules, measurement conflicts, and existing infrastructure tests.

No split, feature filtering, imputation, scaler, or model was fitted.
