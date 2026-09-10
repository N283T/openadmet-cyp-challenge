CREATE TABLE pretraining_sets (
    run_id text PRIMARY KEY REFERENCES runs(id),
    dataset_id text NOT NULL REFERENCES molecular_datasets(run_id),
    split_run_id text NOT NULL REFERENCES split_sets(run_id),
    outer_fold integer NOT NULL CHECK (outer_fold >= 0),
    initialization text NOT NULL,
    summary jsonb NOT NULL
);
CREATE TABLE pretraining_membership (
    run_id text NOT NULL REFERENCES pretraining_sets(run_id),
    compound_id text NOT NULL REFERENCES compounds(id),
    role text NOT NULL CHECK (role IN ('inner_train', 'inner_validation', 'outer_holdout')),
    PRIMARY KEY (run_id, compound_id)
);
CREATE TABLE auxiliary_predictions (
    run_id text NOT NULL REFERENCES pretraining_sets(run_id),
    compound_id text NOT NULL REFERENCES compounds(id),
    endpoint text NOT NULL,
    fold integer NOT NULL CHECK (fold >= 0),
    y_true double precision NOT NULL,
    prediction double precision NOT NULL,
    PRIMARY KEY (run_id, compound_id, endpoint)
);
