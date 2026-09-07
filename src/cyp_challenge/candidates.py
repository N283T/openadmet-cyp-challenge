"""Freeze candidate files and prepare separate submission copies."""

import shutil
from pathlib import Path

from psycopg.types.json import Jsonb

from .database import connect
from .storage import allocate, digest, now, write_json
from .tracks import validate


def candidate(args):
    with connect() as conn:
        recorded_hashes = set()
        for run_id in args.run:
            row = conn.execute(
                "SELECT track,status,manifest FROM runs WHERE id=%s", (run_id,)
            ).fetchone()
            if not row or row[1] != "completed" or row[0] not in (args.track, "shared"):
                raise ValueError(f"Run must be completed and track-compatible: {run_id}")
            recorded_hashes.update(
                item["sha256"]
                for item in row[2].get("files", [])
                if item["path"].startswith("outputs/")
            )
        if digest(args.file) not in recorded_hashes:
            raise ValueError("Prediction file is not a recorded output of the cited runs")
        identifier, path = allocate("candidates", args.track)
        target = path / ("predictions" + args.file.suffix)
        shutil.copy2(args.file, target)
        if digest(target) not in recorded_hashes:
            raise ValueError("Prediction file changed while copying")
        reference = path / ("reference" + args.reference.suffix)
        shutil.copy2(args.reference, reference)
        report = validate(target, args.track, reference)
        manifest = {
            "id": identifier,
            "track": args.track,
            "created_at": now(),
            "runs": args.run,
            "notes": args.notes,
            "file": str(target),
            "sha256": digest(target),
            "reference_sha256": digest(reference),
            "validation": report,
        }
        write_json(path / "manifest.json", manifest)
        conn.execute(
            "INSERT INTO candidates(id,track,manifest_path,manifest) VALUES (%s,%s,%s,%s)",
            (identifier, args.track, str(path / "manifest.json"), Jsonb(manifest)),
        )
        for run_id in set(args.run):
            conn.execute("INSERT INTO candidate_runs VALUES (%s,%s)", (identifier, run_id))
    target.chmod(0o444)
    print(identifier, path)


def prepare(identifier):
    with connect() as conn:
        row = conn.execute(
            "SELECT track,manifest FROM candidates WHERE id=%s", (identifier,)
        ).fetchone()
        if not row:
            raise ValueError("Unknown candidate")
        track, source = row
        if digest(source["file"]) != source["sha256"]:
            raise ValueError("Candidate file was modified")
        submission_id, path = allocate("submissions", track)
        file = path / Path(source["file"]).name
        shutil.copy2(source["file"], file)
        if digest(file) != source["sha256"]:
            raise ValueError("Candidate changed while copying")
        manifest = {
            "id": submission_id,
            "candidate_id": identifier,
            "track": track,
            "status": "prepared",
            "created_at": now(),
            "file": str(file),
            "sha256": digest(file),
        }
        write_json(path / "manifest.json", manifest)
        conn.execute(
            "INSERT INTO submissions(id,candidate_id,track,status,manifest_path,manifest) "
            "VALUES (%s,%s,%s,'prepared',%s,%s)",
            (submission_id, identifier, track, str(path / "manifest.json"), Jsonb(manifest)),
        )
    print(submission_id, path)


def event(args):
    with connect() as conn:
        row = conn.execute(
            "SELECT manifest_path FROM submissions WHERE id=%s", (args.id,)
        ).fetchone()
        if not row:
            raise ValueError("Unknown submission")
        from uuid import uuid4

        path = Path(row[0]).parent / "events" / uuid4().hex
        path.mkdir(parents=True)
        evidence = path / args.evidence.name
        shutil.copy2(args.evidence, evidence)
        record = {
            "status": args.status,
            "notes": args.notes,
            "created_at": now(),
            "evidence_path": str(evidence),
            "sha256": digest(evidence),
        }
        write_json(path / "event.json", record)
        conn.execute(
            "INSERT INTO submission_events(submission_id,status,evidence_path,"
            "evidence_sha256,notes) VALUES (%s,%s,%s,%s,%s)",
            (args.id, args.status, str(evidence), record["sha256"], args.notes),
        )
        conn.execute("UPDATE submissions SET status=%s WHERE id=%s", (args.status, args.id))
    print(path)
