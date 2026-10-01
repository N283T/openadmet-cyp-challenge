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

from . import storage
from .database import connect
from .storage import ROOT, allocate, digest, inventory, now, snapshot_source, write_json


def resolve_inputs(specs, replay=False):
    """Resolve inputs before allocating a run and assign stable snapshot names."""
    resolved = []
    names = set()
    for index, spec in enumerate(specs):
        raw = os.fspath(spec)
        alias = None
        source_text = raw
        if not replay and "=" in raw:
            alias, source_text = raw.split("=", 1)
            if not alias or Path(alias).name != alias or alias in {".", ".."}:
                raise ValueError(f"Invalid input alias: {alias!r}")
        source = Path(source_text).resolve(strict=True)
        name = source.name if replay else alias or f"{index:03d}-{source.name}"
        if name in names:
            raise ValueError(f"Duplicate input snapshot name: {name}")
        names.add(name)
        resolved.append((source, name))
    return resolved


def snapshot_input(source, target, mode="copy"):
    """Link retained immutable assets explicitly; never chmod the shared inode."""
    if mode == "hardlink":
        os.link(source, target)
    elif mode == "copy":
        shutil.copy2(source, target)
        target.chmod(0o444)
    else:
        raise ValueError(f"Unknown input storage mode: {mode}")


def run(args):
    command = args.command
    if command and command[0] == "--":
        command = command[1:]
    if not command:
        raise ValueError("An explicit command is required")
    replay = bool(getattr(args, "replay_of", None))
    resolved_inputs = resolve_inputs(args.input, replay=replay)
    config = args.config.resolve(strict=True) if args.config else None
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
        for source, name in resolved_inputs:
            target = path / "inputs" / name
            mode = getattr(args, "input_mode", "copy")
            snapshot_input(source, target, mode)
            manifest["inputs"].append(
                {
                    "original": str(source),
                    "snapshot": str(target),
                    "sha256": digest(target),
                    "storage": mode,
                }
            )
        if config:
            target = path / ("config" + config.suffix)
            shutil.copy2(config, target)
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
        hashes = {entry["path"]: entry["sha256"] for entry in result["files"]}
        for item in manifest["inputs"]:
            key = str(Path(item["snapshot"]).relative_to(path))
            if hashes.get(key) != item["sha256"]:
                result.update(
                    status="failed", exit_code=1, error=f"Input changed during run: {key}"
                )
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
            input_mode="hardlink"
            if all(v.get("storage") == "hardlink" for v in manifest["inputs"])
            else "copy",
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


def _directory_bytes(path):
    return sum(entry.stat().st_size for entry in path.rglob("*") if entry.is_file())


def _read_json(path):
    try:
        return json.loads(path.read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        return {}


def _console_tail(path, limit=4096):
    try:
        with path.open("rb") as stream:
            stream.seek(0, os.SEEK_END)
            size = stream.tell()
            stream.seek(max(0, size - limit))
            return stream.read().decode(errors="replace")
    except FileNotFoundError:
        return ""


def prune(args):
    """Remove failed run payloads while retaining a compact failure tombstone."""
    with connect() as conn:
        rows = conn.execute(
            "SELECT id,track,purpose,status,issue,manifest_path,created_at,finished_at "
            "FROM runs WHERE status=%s ORDER BY created_at",
            (args.status,),
        ).fetchall()
        note_rows = (
            conn.execute(
                "SELECT run_id,decision,notes,issue,created_at FROM run_notes "
                "WHERE run_id = ANY(%s) ORDER BY created_at",
                ([row[0] for row in rows],),
            ).fetchall()
            if rows
            else []
        )

    notes_by_run = {}
    for run_id, decision, notes, issue, created_at in note_rows:
        notes_by_run.setdefault(run_id, []).append(
            {
                "decision": decision,
                "notes": notes,
                "issue": issue,
                "created_at": created_at.isoformat(),
            }
        )

    planned = []
    runs_root = (storage.ROOT / "artifacts" / "runs").resolve()
    for identifier, track, purpose, status, issue, manifest_path, created_at, finished_at in rows:
        manifest_file = Path(manifest_path).resolve()
        tombstone_file = (
            storage.ROOT / "artifacts" / "pruned-runs" / track / f"{identifier}.json"
        ).resolve()
        if manifest_file == tombstone_file:
            continue
        path = manifest_file.parent
        expected = (runs_root / track / identifier).resolve()
        if path != expected or path.parent.parent != runs_root:
            raise ValueError(f"Refusing unexpected run path for {identifier}: {path}")
        if path.is_dir():
            planned.append(
                (identifier, track, purpose, status, issue, path, created_at, finished_at)
            )

    sizes = {identifier: _directory_bytes(path) for identifier, _, _, _, _, path, _, _ in planned}
    total = sum(sizes.values())
    mode = "PRUNE" if args.execute else "DRY-RUN"
    for identifier, track, purpose, _, _, _, _, _ in planned:
        print(f"{mode}\t{sizes[identifier]}\t{track}\t{identifier}\t{purpose}")
    print(f"{mode} runs={len(planned)} bytes={total}")
    if not args.execute:
        return

    pruned_at = now()
    for identifier, track, purpose, status, issue, path, created_at, finished_at in planned:
        intent = _read_json(path / "intent.json")
        result = _read_json(path / "result.json")
        tombstone_path = storage.ROOT / "artifacts" / "pruned-runs" / track / f"{identifier}.json"
        tombstone_path.parent.mkdir(parents=True, exist_ok=True)
        tombstone = {
            "id": identifier,
            "track": track,
            "purpose": purpose,
            "status": status,
            "issue": issue,
            "created_at": created_at.isoformat(),
            "finished_at": finished_at.isoformat() if finished_at else result.get("finished_at"),
            "command": intent.get("command"),
            "failure": {
                "exit_code": result.get("exit_code"),
                "error": result.get("error"),
                "console_tail": _console_tail(path / "console.log"),
            },
            "notes": notes_by_run.get(identifier, []),
            "pruned_at": pruned_at,
            "prune_reason": "Failed run payload removed by retention policy",
            "reclaimed_bytes": sizes[identifier],
            "original_path": str(path),
            "tombstone_path": str(tombstone_path),
        }
        write_json(tombstone_path, tombstone)
        with connect() as conn:
            conn.execute(
                "UPDATE runs SET manifest_path=%s,manifest=%s WHERE id=%s AND status=%s",
                (str(tombstone_path), Jsonb(tombstone), identifier, status),
            )
            conn.execute(
                "UPDATE run_notes SET record_path=%s WHERE run_id=%s",
                (str(tombstone_path), identifier),
            )
        shutil.rmtree(path)

    print(f"PRUNED runs={len(planned)} bytes={total}")


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
