CREATE TABLE submission_attempts (
 id text PRIMARY KEY,
 submission_id text NOT NULL REFERENCES submissions(id),
 track text NOT NULL,
 file_sha256 text NOT NULL,
 status text NOT NULL CHECK (status IN ('unknown','sent','rejected')),
 manifest_path text NOT NULL,
 created_at timestamptz NOT NULL DEFAULT now()
);
CREATE INDEX submission_attempt_payload ON submission_attempts(track,file_sha256);
