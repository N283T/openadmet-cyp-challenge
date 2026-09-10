CREATE TABLE readout_sets (
    run_id text PRIMARY KEY REFERENCES runs(id),
    dataset_id text NOT NULL REFERENCES molecular_datasets(run_id),
    encoder_run_id text NOT NULL REFERENCES pretraining_sets(run_id),
    split_run_id text NOT NULL REFERENCES split_sets(run_id),
    fold integer NOT NULL CHECK (fold >= 0),
    summary jsonb NOT NULL
);
CREATE TABLE readout_predictions (
    run_id text NOT NULL REFERENCES readout_sets(run_id),
    compound_id text NOT NULL REFERENCES compounds(id),
    endpoint text NOT NULL,
    fold integer NOT NULL CHECK (fold >= 0),
    y_true double precision NOT NULL,
    prediction double precision NOT NULL,
    PRIMARY KEY (run_id, compound_id, endpoint)
);
