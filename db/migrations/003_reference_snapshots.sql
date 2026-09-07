CREATE TABLE reference_snapshots (
    id text PRIMARY KEY,
    url text NOT NULL,
    path text NOT NULL UNIQUE,
    metadata jsonb NOT NULL,
    created_at timestamptz NOT NULL DEFAULT now()
);
