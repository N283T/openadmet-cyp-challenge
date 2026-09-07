CREATE TABLE run_notes (
    id text PRIMARY KEY,
    run_id text NOT NULL REFERENCES runs(id),
    decision text NOT NULL CHECK (decision IN ('observation','keep','discard')),
    notes text NOT NULL,
    issue text,
    record_path text NOT NULL,
    created_at timestamptz NOT NULL DEFAULT now()
);
