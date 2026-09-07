"""Durable, append-only artifact directories and provenance."""

import hashlib
import json
import shutil
import subprocess
import uuid
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def now():
    return datetime.now(UTC).isoformat()


def digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def allocate(kind, track=None):
    identifier = datetime.now(UTC).strftime("%Y%m%dT%H%M%S%fZ") + "-" + uuid.uuid4().hex[:8]
    parent = ROOT / "artifacts" / kind
    if track:
        parent /= track
    path = parent / identifier
    path.mkdir(parents=True, exist_ok=False)
    return identifier, path


def write_json(path, value):
    """Never replace an existing record."""
    with Path(path).open("x") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")


def inventory(path):
    return [
        {"path": str(p.relative_to(path)), "sha256": digest(p), "bytes": p.stat().st_size}
        for p in sorted(Path(path).rglob("*"))
        if p.is_file() and not p.is_symlink()
    ]


def snapshot_source(destination):
    destination.mkdir()
    for name in ("src", "configs", "scripts", "db/migrations"):
        source = ROOT / name
        if source.exists():
            shutil.copytree(
                source, destination / name, ignore=shutil.ignore_patterns("__pycache__")
            )
    for name in ("pixi.toml", "pixi.lock", "pyproject.toml", "AGENTS.md"):
        if (ROOT / name).exists():
            shutil.copy2(ROOT / name, destination / name)

    def git(*args):
        result = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True)
        return result.stdout.strip() if result.returncode == 0 else None

    (destination / "working-tree.patch").write_text(git("diff", "HEAD", "--binary") or "")
    return {
        "commit": git("rev-parse", "HEAD"),
        "status": git("status", "--porcelain"),
        "source_files": inventory(destination),
    }
