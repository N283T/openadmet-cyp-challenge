CREATE EXTENSION IF NOT EXISTS rdkit;
CREATE TABLE runs (
    id text PRIMARY KEY,
    track text NOT NULL CHECK (track IN ('shared', 'regression', 'classification')),
    purpose text NOT NULL,
    status text NOT NULL CHECK (status IN ('running','completed','failed','interrupted')),
    issue text,
    manifest_path text NOT NULL UNIQUE,
    manifest jsonb NOT NULL,
    created_at timestamptz NOT NULL DEFAULT now(),
    finished_at timestamptz
);
CREATE TABLE assets (
    id text PRIMARY KEY,
    kind text NOT NULL CHECK (kind IN ('raw','api')),
    manifest_path text NOT NULL UNIQUE,
    manifest jsonb NOT NULL,
    created_at timestamptz NOT NULL DEFAULT now()
);
CREATE TABLE candidates (
    id text PRIMARY KEY,
    track text NOT NULL CHECK (track IN ('regression','classification')),
    manifest_path text NOT NULL UNIQUE,
    manifest jsonb NOT NULL,
    created_at timestamptz NOT NULL DEFAULT now()
);
CREATE TABLE candidate_runs (
    candidate_id text REFERENCES candidates(id),
    run_id text REFERENCES runs(id),
    PRIMARY KEY(candidate_id, run_id)
);
CREATE TABLE submissions (
    id text PRIMARY KEY,
    candidate_id text NOT NULL REFERENCES candidates(id),
    track text NOT NULL CHECK (track IN ('regression','classification')),
    status text NOT NULL CHECK (status IN ('prepared','sent','unknown','accepted','rejected')),
    manifest_path text NOT NULL UNIQUE,
    manifest jsonb NOT NULL,
    created_at timestamptz NOT NULL DEFAULT now()
);
CREATE TABLE submission_events (
    id bigint GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    submission_id text NOT NULL REFERENCES submissions(id),
    status text NOT NULL CHECK (status IN ('sent','unknown','accepted','rejected')),
    evidence_path text NOT NULL,
    evidence_sha256 text NOT NULL,
    notes text NOT NULL,
    created_at timestamptz NOT NULL DEFAULT now()
);
CREATE TABLE metrics (
    run_id text NOT NULL REFERENCES runs(id),
    endpoint text NOT NULL,
    split_id text NOT NULL,
    fold text NOT NULL,
    metric text NOT NULL,
    value double precision NOT NULL,
    PRIMARY KEY (run_id, endpoint, split_id, fold, metric)
);
