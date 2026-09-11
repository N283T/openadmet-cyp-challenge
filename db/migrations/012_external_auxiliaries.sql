CREATE TABLE external_datasets (
    run_id text PRIMARY KEY REFERENCES runs(id),
    chembl_run_id text NOT NULL REFERENCES runs(id),
    pubchem_run_id text NOT NULL REFERENCES runs(id),
    summary jsonb NOT NULL,
    loaded_at timestamptz NOT NULL DEFAULT now()
);
CREATE TABLE external_compounds (
    dataset_id text NOT NULL REFERENCES external_datasets(run_id),
    compound_id text NOT NULL REFERENCES compounds(id),
    connectivity text NOT NULL,
    PRIMARY KEY (dataset_id, compound_id)
);
CREATE TABLE external_targets (
    dataset_id text NOT NULL REFERENCES external_datasets(run_id),
    compound_id text NOT NULL REFERENCES compounds(id),
    endpoint text NOT NULL,
    source text NOT NULL,
    value double precision NOT NULL,
    source_rows integer NOT NULL CHECK (source_rows > 0),
    assays integer NOT NULL CHECK (assays > 0),
    PRIMARY KEY (dataset_id, compound_id, endpoint)
);
CREATE TABLE external_target_sources (
    dataset_id text NOT NULL,
    compound_id text NOT NULL,
    endpoint text NOT NULL,
    source_record_id text NOT NULL,
    raw_value double precision NOT NULL,
    metadata jsonb NOT NULL,
    PRIMARY KEY (dataset_id, endpoint, source_record_id),
    FOREIGN KEY (dataset_id, compound_id, endpoint)
        REFERENCES external_targets(dataset_id, compound_id, endpoint)
);
