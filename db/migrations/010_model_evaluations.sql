CREATE TABLE model_evaluations (
 run_id text PRIMARY KEY REFERENCES runs(id),
 audit_run_id text NOT NULL REFERENCES runs(id),
 dataset_run_id text NOT NULL REFERENCES molecular_datasets(run_id),
 split_run_id text NOT NULL REFERENCES split_sets(run_id),
 variant text NOT NULL,
 fold integer NOT NULL CHECK (fold BETWEEN 0 AND 4),
 summary jsonb NOT NULL
);
CREATE TABLE evaluation_predictions (
 run_id text NOT NULL REFERENCES model_evaluations(run_id),
 compound_id text NOT NULL REFERENCES compounds(id),
 endpoint text NOT NULL,
 fold integer NOT NULL,
 y_true double precision NOT NULL,
 prediction double precision NOT NULL,
 PRIMARY KEY (run_id,compound_id,endpoint)
);
