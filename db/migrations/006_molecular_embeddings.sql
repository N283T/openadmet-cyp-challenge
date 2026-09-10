CREATE TABLE embedding_sets (
    run_id text PRIMARY KEY REFERENCES runs(id),
    dataset_id text NOT NULL REFERENCES molecular_datasets(run_id),
    family text NOT NULL,
    dimension integer NOT NULL CHECK (dimension > 0),
    summary jsonb NOT NULL,
    UNIQUE (run_id, dimension)
);
CREATE TABLE compound_embeddings (
    feature_set_id text NOT NULL,
    compound_id text NOT NULL REFERENCES compounds(id),
    dimension integer NOT NULL,
    embedding real[] NOT NULL,
    PRIMARY KEY (feature_set_id, compound_id),
    FOREIGN KEY (feature_set_id, dimension) REFERENCES embedding_sets(run_id, dimension),
    CHECK (array_ndims(embedding) = 1 AND cardinality(embedding) = dimension),
    CHECK (array_position(embedding, NULL) IS NULL)
);
