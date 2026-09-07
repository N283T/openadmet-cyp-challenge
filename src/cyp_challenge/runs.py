"""Record exploratory and production commands before execution."""

import json
import os
import platform
import shutil
import signal
import subprocess
import sys
from pathlib import Path

from psycopg.types.json import Jsonb

from .database import connect
from .storage import ROOT, allocate, digest, inventory, now, snapshot_source, write_json


def run(args):
    command = args.command
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        raise ValueError("An explicit command is required")
    with connect():
        pass
    identifier, path = allocate("runs", args.track)
    manifest = {
        "id": identifier,
        "track": args.track,
        "purpose": args.purpose,
        "issue": args.issue,
        "created_at": now(),
        "command": command,
        "cwd": str(path / "source"),
        "status": "running",
        "inputs": [],
        "python": sys.version,
        "platform": platform.platform(),
        "replay_of": getattr(args, "replay_of", None),
    }
    write_json(path / "intent.json", manifest)
    with connect() as conn:
        conn.execute(
            "INSERT INTO runs (id,track,purpose,status,issue,manifest_path,manifest) "
            "VALUES (%s,%s,%s,'running',%s,%s,%s)",
            (
                identifier,
                args.track,
                args.purpose,
                args.issue,
                str(path / "manifest.json"),
                Jsonb(manifest),
            ),
        )
    print(f"run_id={identifier}\n{path}", flush=True)
    result = {"id": identifier, "status": "failed", "exit_code": None}
    process = None
    try:
        for name in ("outputs", "scratch", "inputs"):
            (path / name).mkdir()
        for i, source in enumerate(args.input):
            source = source.resolve(strict=True)
            name = source.name if getattr(args, "replay_of", None) else f"{i:03d}-{source.name}"
            target = path / "inputs" / name
            shutil.copy2(source, target)
            target.chmod(0o444)
            manifest["inputs"].append(
                {"original": str(source), "snapshot": str(target), "sha256": digest(target)}
            )
        if args.config:
            target = path / ("config" + args.config.suffix)
            shutil.copy2(args.config, target)
            manifest["config"] = {"path": str(target), "sha256": digest(target)}
        if getattr(args, "source_snapshot", None):
            shutil.copytree(args.source_snapshot, path / "source")
            manifest["git"] = {
                "replayed_source": str(args.source_snapshot),
                "source_files": inventory(path / "source"),
            }
        else:
            manifest["git"] = snapshot_source(path / "source")
        with (path / "environment.json").open("x") as stream:
            subprocess.run(["pixi", "list", "--json"], cwd=ROOT, stdout=stream, check=True)
        write_json(path / "manifest.json", manifest)
        with connect() as conn:
            conn.execute("UPDATE runs SET manifest=%s WHERE id=%s", (Jsonb(manifest), identifier))
        env = os.environ.copy()
        env.update(
            CYP_RUN_DIR=str(path),
            CYP_OUTPUT_DIR=str(path / "outputs"),
            CYP_INPUT_DIR=str(path / "inputs"),
            TMPDIR=str(path / "scratch"),
            PYTHONPATH=str(path / "source/src"),
        )
        if args.config:
            env["CYP_CONFIG_PATH"] = manifest["config"]["path"]
        with (path / "console.log").open("xb") as log:
            process = subprocess.Popen(
                command,
                cwd=path / "source",
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                start_new_session=True,
            )
            write_json(path / "process.json", {"pid": process.pid, "host": platform.node()})
            code = process.wait()
        result.update(status="completed" if code == 0 else "failed", exit_code=code)
    except KeyboardInterrupt:
        if process:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
        result.update(status="interrupted", exit_code=130)
    except Exception as exc:
        result["error"] = repr(exc)
        print(f"Run failed: {exc}", file=sys.stderr)
    finally:
        result.update(finished_at=now(), files=inventory(path))
        write_json(path / "result.json", result)
        with connect() as conn:
            conn.execute(
                "UPDATE runs SET status=%s,finished_at=%s,manifest=manifest || %s WHERE id=%s",
                (result["status"], result["finished_at"], Jsonb(result), identifier),
            )
    return result["exit_code"] if result["exit_code"] is not None else 1


def replay(identifier):
    """Replay verified snapshots with the current identical lock; never mutate the old run."""
    from argparse import Namespace

    with connect() as conn:
        row = conn.execute("SELECT manifest FROM runs WHERE id=%s", (identifier,)).fetchone()
    if not row or row[0]["status"] != "completed":
        raise ValueError("Replay requires a completed run")
    manifest = row[0]
    path = Path(manifest["cwd"]).parent
    for file in manifest["files"]:
        if digest(path / file["path"]) != file["sha256"]:
            raise ValueError(f"Snapshot changed: {file['path']}")
    if digest(path / "source/pixi.lock") != digest(ROOT / "pixi.lock"):
        raise ValueError(
            "Lock differs; restore the recorded environment in a separate checkout first"
        )
    return run(
        Namespace(
            track=manifest["track"],
            purpose=f"Replay of {identifier}",
            issue=manifest["issue"],
            command=manifest["command"],
            input=[Path(v["snapshot"]) for v in manifest["inputs"]],
            config=Path(manifest["config"]["path"]) if "config" in manifest else None,
            replay_of=identifier,
            source_snapshot=path / "source",
        )
    )


def reconcile():
    """Restore completion records after interruption between filesystem and DB writes."""
    with connect() as conn:
        for identifier, manifest_path in conn.execute(
            "SELECT id,manifest_path FROM runs WHERE status='running'"
        ).fetchall():
            result_path = Path(manifest_path).with_name("result.json")
            if result_path.exists():
                result = json.loads(result_path.read_text())
                conn.execute(
                    "UPDATE runs SET status=%s,finished_at=%s,manifest=manifest || %s WHERE id=%s",
                    (result["status"], result["finished_at"], Jsonb(result), identifier),
                )
                print(f"Reconciled {identifier}")
            else:
                print(f"Needs inspection (possibly still running): {identifier} {manifest_path}")


def note(args):
    """Preserve reasons for keeping or discarding even a successful throwaway probe."""
    from uuid import uuid4

    with connect() as conn:
        row = conn.execute("SELECT manifest_path FROM runs WHERE id=%s", (args.id,)).fetchone()
        if not row:
            raise ValueError("Unknown run")
        identifier = uuid4().hex
        path = Path(row[0]).parent / "notes" / identifier
        path.mkdir(parents=True)
        record = {
            "id": identifier,
            "run_id": args.id,
            "decision": args.decision,
            "notes": args.notes,
            "issue": args.issue,
            "created_at": now(),
        }
        write_json(path / "note.json", record)
        conn.execute(
            "INSERT INTO run_notes(id,run_id,decision,notes,issue,record_path) "
            "VALUES (%s,%s,%s,%s,%s,%s)",
            (identifier, args.id, args.decision, args.notes, args.issue, str(path / "note.json")),
        )
    print(path)
