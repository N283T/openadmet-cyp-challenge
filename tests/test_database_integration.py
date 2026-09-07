"""Opt-in integration tests using an isolated schema in the local challenge DB."""

import json
import os
import shutil
from argparse import Namespace
from pathlib import Path
from uuid import uuid4

import polars as pl
import pytest
from psycopg import sql

from cyp_challenge import candidates, database, runs, storage
from cyp_challenge.tracks import ENDPOINTS

pytestmark = pytest.mark.skipif(os.environ.get("CYP_TEST_DB") != "1", reason="Set CYP_TEST_DB=1")


@pytest.fixture
def registry(monkeypatch, tmp_path):
    schema = "test_" + uuid4().hex
    original = database.connect
    with original() as conn:
        conn.execute(sql.SQL("CREATE SCHEMA {}").format(sql.Identifier(schema)))

    def connect():
        connection = original()
        connection.execute(sql.SQL("SET search_path TO {}, public").format(sql.Identifier(schema)))
        return connection

    try:
        with connect() as conn:
            for migration in sorted((storage.ROOT / "db/migrations").glob("*.sql")):
                conn.execute(
                    migration.read_text().replace("CREATE EXTENSION IF NOT EXISTS rdkit;", "")
                )
        for module in (runs, candidates):
            monkeypatch.setattr(module, "connect", connect)
        shutil.copy2(storage.ROOT / "pixi.lock", tmp_path / "pixi.lock")
        monkeypatch.setattr(storage, "ROOT", tmp_path)
        yield connect, tmp_path
    finally:
        with original() as conn:
            conn.execute(sql.SQL("DROP SCHEMA {} CASCADE").format(sql.Identifier(schema)))


def launch(command, **kwargs):
    return runs.run(
        Namespace(
            track="shared",
            purpose="Isolated infrastructure test",
            issue=None,
            config=None,
            input=[],
            command=command,
            **kwargs,
        )
    )


def test_failed_command_and_missing_input_remain_recorded(registry):
    connect, root = registry
    assert launch(["python", "-c", "raise SystemExit(7)"]) == 7
    args = Namespace(
        track="shared",
        purpose="Missing input test",
        issue=None,
        config=None,
        input=[root / "does-not-exist"],
        command=["python", "-c", "pass"],
    )
    assert runs.run(args) == 1
    with connect() as conn:
        rows = conn.execute("SELECT status,manifest_path FROM runs").fetchall()
    assert len(rows) == 2
    for status, path in rows:
        assert status == "failed"
        assert Path(path).with_name("intent.json").exists()
        assert json.loads(Path(path).with_name("result.json").read_text())["status"] == "failed"


def test_replay_and_tampering(registry):
    connect, root = registry
    command = [
        "python",
        "-c",
        "import os; from pathlib import Path; "
        "Path(os.environ['CYP_OUTPUT_DIR'], 'value.txt').write_text('stable')",
    ]
    assert launch(command) == 0
    with connect() as conn:
        identifier, manifest = conn.execute("SELECT id,manifest FROM runs").fetchone()
    assert runs.replay(identifier) == 0
    results = list((root / "artifacts/runs/shared").glob("*/outputs/value.txt"))
    assert len(results) == 2 and len({storage.digest(p) for p in results}) == 1
    original_output = Path(manifest["cwd"]).parent / "outputs/value.txt"
    original_output.write_text("changed")
    with pytest.raises(ValueError, match="Snapshot changed"):
        runs.replay(identifier)


def test_candidate_and_submission_separate_and_verified(registry):
    connect, root = registry
    reference = root / "reference.csv"
    frame = pl.DataFrame({"SMILES": ["C"] * 750, "Molecule_Name": [f"m{i}" for i in range(750)]})
    frame.write_csv(reference)
    for endpoint in ENDPOINTS["regression"]:
        frame = frame.with_columns(pl.Series(endpoint, [float(i) for i in range(750)]))
    file = root / "prediction.parquet"
    frame.write_parquet(file)
    assert (
        launch(
            [
                "python",
                "-c",
                "import shutil,os; from pathlib import Path; "
                f"shutil.copy2({str(file)!r}, "
                "Path(os.environ['CYP_OUTPUT_DIR'])/'predictions.parquet')",
            ]
        )
        == 0
    )
    with connect() as conn:
        run_id = conn.execute("SELECT id FROM runs").fetchone()[0]
    args = Namespace(
        track="regression",
        file=file,
        reference=reference,
        run=[run_id],
        notes="Isolated test; never upload",
    )
    candidates.candidate(args)
    with connect() as conn:
        candidate_id, candidate = conn.execute("SELECT id,manifest FROM candidates").fetchone()
    candidates.prepare(candidate_id)
    with connect() as conn:
        submission = conn.execute("SELECT manifest FROM submissions").fetchone()[0]
    assert submission["status"] == "prepared"
    assert submission["file"] != candidate["file"] != str(file)
    assert storage.digest(submission["file"]) == storage.digest(candidate["file"])
    file.write_bytes(b"original overwritten")
    assert storage.digest(candidate["file"]) == candidate["sha256"]
    frozen = Path(candidate["file"])
    frozen.chmod(0o644)
    frozen.write_bytes(b"candidate corrupted")
    with pytest.raises(ValueError, match="modified"):
        candidates.prepare(candidate_id)


def test_discard_reason_is_append_only(registry):
    connect, root = registry
    assert launch(["python", "-c", "pass"]) == 0
    with connect() as conn:
        identifier = conn.execute("SELECT id FROM runs").fetchone()[0]
    runs.note(Namespace(id=identifier, decision="discard", notes="No gain", issue=None))
    runs.note(Namespace(id=identifier, decision="observation", notes="Revisit later", issue=None))
    with connect() as conn:
        rows = conn.execute("SELECT decision,record_path FROM run_notes").fetchall()
    assert len(rows) == 2
    assert {decision for decision, _ in rows} == {"discard", "observation"}
    assert all(Path(path).is_file() for _, path in rows)
