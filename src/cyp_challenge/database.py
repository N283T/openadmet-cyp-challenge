"""Local PostgreSQL lifecycle and versioned schema."""

import subprocess

import psycopg
from psycopg.types.json import Jsonb

from .storage import ROOT, allocate, digest, inventory, now, write_json

RUNTIME = ROOT / ".runtime"
PGDATA = RUNTIME / "postgres"
SOCKET = RUNTIME / "socket"
PORT = "5434"


def connect():
    return psycopg.connect(host=str(SOCKET), port=PORT, dbname="cyp_challenge")


def lifecycle(action):
    RUNTIME.mkdir(exist_ok=True)
    SOCKET.mkdir(exist_ok=True)
    if action == "init":
        if not (PGDATA / "PG_VERSION").exists():
            subprocess.run(
                [
                    "initdb",
                    "-D",
                    str(PGDATA),
                    "--auth-local=trust",
                    "--auth-host=scram-sha-256",
                    "--encoding=UTF8",
                    "--locale=C.UTF-8",
                ],
                check=True,
            )
        return
    if action == "start":
        status = subprocess.run(["pg_ctl", "-D", str(PGDATA), "status"], capture_output=True)
        if status.returncode:
            subprocess.run(
                [
                    "pg_ctl",
                    "-D",
                    str(PGDATA),
                    "-l",
                    str(RUNTIME / "postgres.log"),
                    "-o",
                    f"-k {SOCKET} -p {PORT} -h ''",
                    "-w",
                    "start",
                ],
                check=True,
            )
        with psycopg.connect(
            host=str(SOCKET), port=PORT, dbname="postgres", autocommit=True
        ) as conn:
            if not conn.execute(
                "SELECT 1 FROM pg_database WHERE datname='cyp_challenge'"
            ).fetchone():
                conn.execute("CREATE DATABASE cyp_challenge")
        return
    if action == "stop":
        subprocess.run(["pg_ctl", "-D", str(PGDATA), "-m", "fast", "-w", "stop"], check=True)
    elif action == "migrate":
        with connect() as conn:
            conn.execute("SELECT pg_advisory_xact_lock(847291)")
            conn.execute(
                "CREATE TABLE IF NOT EXISTS schema_migrations "
                "(name text PRIMARY KEY, sha256 text NOT NULL, "
                "applied_at timestamptz DEFAULT now())"
            )
            for file in sorted((ROOT / "db/migrations").glob("*.sql")):
                existing = conn.execute(
                    "SELECT sha256 FROM schema_migrations WHERE name=%s", (file.name,)
                ).fetchone()
                if existing:
                    if existing[0] != digest(file):
                        raise ValueError(f"Applied migration changed: {file}")
                else:
                    conn.execute(file.read_text())
                    conn.execute(
                        "INSERT INTO schema_migrations(name,sha256) VALUES (%s,%s)",
                        (file.name, digest(file)),
                    )
    elif action == "backup":
        identifier, path = allocate("backups")
        subprocess.run(
            [
                "pg_dump",
                "-h",
                str(SOCKET),
                "-p",
                PORT,
                "-Fc",
                "-f",
                str(path / "database.dump"),
                "cyp_challenge",
            ],
            check=True,
        )
        write_json(
            path / "manifest.json",
            {"id": identifier, "created_at": now(), "files": inventory(path)},
        )
        print(path)


def register_asset(identifier, kind, path, manifest):
    with connect() as conn:
        conn.execute(
            "INSERT INTO assets VALUES (%s,%s,%s,%s,DEFAULT)",
            (identifier, kind, str(path / "manifest.json"), Jsonb(manifest)),
        )
