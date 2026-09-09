CREATE TABLE split_sets (
    run_id text PRIMARY KEY REFERENCES runs(id),
    dataset_id text NOT NULL REFERENCES molecular_datasets(run_id),
    summary jsonb NOT NULL
);
CREATE TABLE split_assignments (
    split_run_id text NOT NULL REFERENCES split_sets(run_id),
    strategy text NOT NULL,
    compound_id text NOT NULL REFERENCES compounds(id),
    fold integer NOT NULL CHECK (fold >= 0),
    group_label text NOT NULL,
    PRIMARY KEY (split_run_id, strategy, compound_id)
);
CREATE TABLE baseline_sets (
    run_id text PRIMARY KEY REFERENCES runs(id),
    split_run_id text NOT NULL REFERENCES split_sets(run_id),
    feature_run_id text NOT NULL REFERENCES feature_sets(run_id),
    summary jsonb NOT NULL
);
CREATE TABLE baseline_predictions (
    run_id text NOT NULL REFERENCES baseline_sets(run_id),
    strategy text NOT NULL,
    features text NOT NULL,
    compound_id text NOT NULL REFERENCES compounds(id),
    endpoint text NOT NULL,
    fold integer NOT NULL CHECK (fold >= 0),
    y_true double precision NOT NULL,
    prediction double precision NOT NULL,
    PRIMARY KEY (run_id, strategy, features, compound_id, endpoint)
);
