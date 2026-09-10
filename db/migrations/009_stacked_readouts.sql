CREATE TABLE stacked_readouts (
    run_id text PRIMARY KEY REFERENCES runs(id),
    split_run_id text NOT NULL REFERENCES split_sets(run_id),
    fold integer NOT NULL CHECK (fold BETWEEN 0 AND 4),
    summary jsonb NOT NULL
);
CREATE TABLE stacked_predictions (
    run_id text NOT NULL REFERENCES stacked_readouts(run_id),
    compound_id text NOT NULL REFERENCES compounds(id),
    endpoint text NOT NULL,
    fold integer NOT NULL,
    y_true double precision NOT NULL,
    prediction double precision NOT NULL,
    PRIMARY KEY (run_id, compound_id, endpoint)
);
CREATE TABLE stacked_auxiliary_parents (
    run_id text NOT NULL REFERENCES stacked_readouts(run_id),
    auxiliary_run_id text NOT NULL REFERENCES runs(id),
    PRIMARY KEY (run_id, auxiliary_run_id)
);
