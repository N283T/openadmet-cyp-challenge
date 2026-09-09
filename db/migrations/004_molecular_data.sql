CREATE TABLE molecular_datasets (
    run_id text PRIMARY KEY REFERENCES runs(id),
    summary jsonb NOT NULL,
    loaded_at timestamptz NOT NULL DEFAULT now()
);
CREATE TABLE compounds (
    id text PRIMARY KEY,
    standardized_smiles text NOT NULL UNIQUE,
    mol mol GENERATED ALWAYS AS (mol_from_smiles(standardized_smiles::cstring)) STORED
);
CREATE INDEX compounds_mol_gist ON compounds USING gist(mol);
CREATE TABLE dataset_compounds (
    dataset_id text REFERENCES molecular_datasets(run_id),
    compound_id text REFERENCES compounds(id),
    PRIMARY KEY(dataset_id, compound_id)
);
CREATE TABLE source_records (
    dataset_id text REFERENCES molecular_datasets(run_id),
    record_id text NOT NULL,
    source_file text NOT NULL,
    role text NOT NULL CHECK(role IN ('train','test','aux')),
    molecule_name text NOT NULL,
    raw_smiles text NOT NULL,
    compound_id text REFERENCES compounds(id),
    payload jsonb NOT NULL,
    PRIMARY KEY(dataset_id, record_id)
);
CREATE INDEX source_records_compound ON source_records(compound_id);
CREATE TABLE record_origins (
    dataset_id text NOT NULL,
    source_file text NOT NULL,
    source_row integer NOT NULL,
    record_id text NOT NULL,
    PRIMARY KEY(dataset_id, source_file, source_row),
    FOREIGN KEY(dataset_id, record_id) REFERENCES source_records(dataset_id, record_id)
);
CREATE TABLE measurements (
    dataset_id text REFERENCES molecular_datasets(run_id),
    measurement_id text NOT NULL,
    compound_id text NOT NULL REFERENCES compounds(id),
    endpoint text NOT NULL,
    value double precision NOT NULL,
    conf_low double precision,
    conf_high double precision,
    std double precision,
    PRIMARY KEY(dataset_id, measurement_id)
);
CREATE INDEX measurements_endpoint_compound ON measurements(endpoint, compound_id);
CREATE TABLE measurement_sources (
    dataset_id text NOT NULL,
    measurement_id text NOT NULL,
    record_id text NOT NULL,
    PRIMARY KEY(dataset_id, measurement_id, record_id),
    FOREIGN KEY(dataset_id, measurement_id) REFERENCES measurements(dataset_id, measurement_id),
    FOREIGN KEY(dataset_id, record_id) REFERENCES source_records(dataset_id, record_id)
);
CREATE TABLE feature_sets (
    run_id text PRIMARY KEY REFERENCES runs(id),
    dataset_id text NOT NULL REFERENCES molecular_datasets(run_id),
    catalog jsonb NOT NULL,
    summary jsonb NOT NULL,
    loaded_at timestamptz NOT NULL DEFAULT now()
);
CREATE TABLE compound_descriptors_full (
    feature_set_id text REFERENCES feature_sets(run_id),
    compound_id text REFERENCES compounds(id),
    descriptors jsonb NOT NULL,
    PRIMARY KEY(feature_set_id, compound_id)
);
CREATE TABLE compound_mordred (
    feature_set_id text REFERENCES feature_sets(run_id),
    compound_id text REFERENCES compounds(id),
    descriptors jsonb NOT NULL,
    PRIMARY KEY(feature_set_id, compound_id)
);
CREATE TABLE descriptor_errors (
    feature_set_id text REFERENCES feature_sets(run_id),
    compound_id text REFERENCES compounds(id),
    errors jsonb NOT NULL,
    PRIMARY KEY(feature_set_id, compound_id)
);
CREATE VIEW single_concentration AS
SELECT dataset_id, record_id, compound_id, molecule_name,
       payload->>'enzyme' AS enzyme,
       (payload->>'concentration_M')::double precision AS concentration_m,
       (payload->>'log2fc_estimate')::double precision AS log2fc_estimate,
       (payload->>'log2fc_std_error')::double precision AS log2fc_std_error,
       (payload->>'log2fc_fdr')::double precision AS log2fc_fdr,
       payload
FROM source_records WHERE role='aux';
