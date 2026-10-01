"""Consolidate completed Boltz trunk tensors into one indexed artifact tree."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import tarfile
from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = ROOT / "artifacts"
DESTINATION = ARTIFACTS / "boltz-trunk" / "post-r1"
NAME_PATTERN = re.compile(r"embeddings_record_(\d{6})\.npz")
ARCHIVE_METADATA = {"trunk_manifest.json", "trunk_run_config.json", "colab_chunk_manifest.json"}

TARGETS = {
    "cyp1a2": {
        "expected": 2163,
        "preparation_run": "20260929T034044572251Z-b2e9435f",
        "archives": ARTIFACTS / "boltz-trunk-colab" / "cyp1a2-r1-a100",
        "local_runs": [],
    },
    "cyp2c9": {
        "expected": 2035,
        "preparation_run": "20260929T034046886086Z-b2f6b507",
        "archives": None,
        "local_runs": ["20260929T062353858054Z-c38fada2"],
    },
    "cyp2d6": {
        "expected": 2247,
        "preparation_run": "20260929T034049410496Z-672bc9b5",
        "archives": None,
        "local_runs": ["20260929T110347354126Z-a56c2a81"],
    },
    "cyp3a4": {
        "expected": 4334,
        "preparation_run": "20260929T034052097919Z-b106051b",
        "archives": ARTIFACTS / "boltz-trunk-colab" / "cyp3a4-r1-a100",
        "local_runs": [
            "20260929T201936481035Z-7e76198a",
            "20260929T160603505210Z-7c34d9eb",
        ],
    },
}


def sha256(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def archive_members(path: Path) -> list[tarfile.TarInfo]:
    with tarfile.open(path, "r:") as archive:
        all_members = archive.getmembers()
    for member in all_members:
        if not member.isfile() or not (
            NAME_PATTERN.fullmatch(member.name) or member.name in ARCHIVE_METADATA
        ):
            raise ValueError(f"Unexpected tar member in {path}: {member.name}")
    return [member for member in all_members if NAME_PATTERN.fullmatch(member.name)]


def archive_digest(path: Path) -> str:
    sidecar = path.with_suffix(path.suffix + ".sha256")
    expected, recorded_name = sidecar.read_text().split()
    if recorded_name != path.name:
        raise ValueError(f"Archive sidecar name mismatch: {sidecar}")
    actual = sha256(path)
    if actual != expected:
        raise ValueError(f"Archive checksum mismatch: {path}")
    return actual


def result_hashes(run_id: str) -> dict[str, str]:
    path = ARTIFACTS / "runs" / "shared" / run_id / "result.json"
    result = json.loads(path.read_text())
    if result.get("status") != "completed":
        raise ValueError(f"Source run is not completed: {run_id}")
    return {
        Path(entry["path"]).name: entry["sha256"]
        for entry in result["files"]
        if NAME_PATTERN.fullmatch(Path(entry["path"]).name)
    }


def discover_sources() -> tuple[dict[str, dict[int, dict]], dict[str, list[dict]]]:
    sources: dict[str, dict[int, dict]] = {}
    archives_by_target: dict[str, list[dict]] = {}
    for target, config in TARGETS.items():
        target_sources: dict[int, dict] = {}
        archive_rows = []
        archive_dir = config["archives"]
        if archive_dir and archive_dir.exists():
            for archive in sorted(archive_dir.glob("chunk-*.tar")):
                members = archive_members(archive)
                indices = []
                for member in members:
                    index = int(NAME_PATTERN.fullmatch(member.name).group(1))
                    if index in target_sources:
                        raise ValueError(f"Duplicate {target} record {index}")
                    target_sources[index] = {
                        "kind": "archive",
                        "path": archive,
                        "member": member.name,
                        "bytes": member.size,
                    }
                    indices.append(index)
                archive_rows.append(
                    {
                        "path": archive,
                        "members": members,
                        "first": min(indices),
                        "last": max(indices),
                    }
                )

        for run_id in config["local_runs"]:
            output = ARTIFACTS / "runs" / "shared" / run_id / "outputs"
            hashes = result_hashes(run_id)
            for path in sorted(output.glob("embeddings_record_*.npz")):
                index = int(NAME_PATTERN.fullmatch(path.name).group(1))
                if index in target_sources:
                    raise ValueError(f"Duplicate {target} record {index}")
                target_sources[index] = {
                    "kind": "run",
                    "path": path,
                    "run_id": run_id,
                    "bytes": path.stat().st_size,
                    "sha256": hashes[path.name],
                }

        records_dir = DESTINATION / target / "records"
        existing_manifest = DESTINATION / target / "manifest.parquet"
        existing_hashes = {}
        if existing_manifest.exists():
            existing_hashes = {
                row["record_index"]: row["sha256"]
                for row in pl.read_parquet(existing_manifest).iter_rows(named=True)
            }
        for path in sorted(records_dir.glob("embeddings_record_*.npz")):
            index = int(NAME_PATTERN.fullmatch(path.name).group(1))
            if index in target_sources:
                continue
            if target == "cyp1a2" or (target == "cyp3a4" and index <= 2000):
                start = ((index - 1) // 250) * 250 + 1
                limit = config["expected"] if target == "cyp1a2" else 2000
                end = min(start + 249, limit)
                target_sources[index] = {
                    "kind": "archive",
                    "path": Path(f"chunk-{start:06d}-{end:06d}.tar"),
                    "member": path.name,
                    "bytes": path.stat().st_size,
                    "sha256": existing_hashes.get(index),
                }
            else:
                raise ValueError(f"Canonical file has no recognized source: {target} {index}")

        expected = set(range(1, config["expected"] + 1))
        actual = set(target_sources)
        if actual != expected:
            missing = sorted(expected - actual)[:10]
            extra = sorted(actual - expected)[:10]
            raise ValueError(f"{target} coverage mismatch: missing={missing}, extra={extra}")
        sources[target] = target_sources
        archives_by_target[target] = archive_rows
    return sources, archives_by_target


def preparation_manifest(target: str) -> pl.DataFrame:
    run_id = TARGETS[target]["preparation_run"]
    path = ARTIFACTS / "runs" / "shared" / run_id / "outputs" / "manifest.parquet"
    frame = pl.read_parquet(path).with_row_index("record_index", offset=1)
    if frame.height != TARGETS[target]["expected"]:
        raise ValueError(f"Preparation manifest count mismatch: {target}")
    prefixes = frame.get_column("record_id").str.slice(0, 6).cast(pl.Int64)
    if prefixes.to_list() != list(range(1, frame.height + 1)):
        raise ValueError(f"Preparation manifest order mismatch: {target}")
    return frame


def extract_archive(archive_row: dict, records_dir: Path) -> dict[str, str]:
    archive = archive_row["path"]
    digest = archive_digest(archive)
    with tarfile.open(archive, "r:") as stream:
        for member in archive_row["members"]:
            destination = records_dir / member.name
            if destination.exists():
                if destination.stat().st_size != member.size:
                    raise ValueError(f"Existing extracted file has wrong size: {destination}")
                continue
            partial = destination.with_suffix(".npz.partial")
            source = stream.extractfile(member)
            if source is None:
                raise ValueError(f"Cannot extract {member.name} from {archive}")
            with source, partial.open("xb") as output:
                shutil.copyfileobj(source, output, length=16 * 1024 * 1024)
            if partial.stat().st_size != member.size:
                raise ValueError(f"Extracted size mismatch: {member.name}")
            partial.replace(destination)
    return {
        "archive": archive.name,
        "sha256": digest,
        "first_record": archive_row["first"],
        "last_record": archive_row["last"],
    }


def move_local(source: Path, destination: Path) -> None:
    if source.is_symlink():
        if source.resolve() != destination.resolve():
            raise ValueError(f"Unexpected existing symlink: {source}")
        return
    if destination.exists():
        raise FileExistsError(destination)
    source.replace(destination)
    source.symlink_to(os.path.relpath(destination, source.parent))


def organize(execute: bool) -> None:
    sources, archives_by_target = discover_sources()
    for target, target_sources in sources.items():
        archive_count = sum(row["kind"] == "archive" for row in target_sources.values())
        local_count = len(target_sources) - archive_count
        print(
            f"{target}: records={len(target_sources)} archive={archive_count} local={local_count}"
        )
    if not execute:
        print(f"DRY-RUN destination={DESTINATION}")
        return

    DESTINATION.mkdir(parents=True, exist_ok=True)
    all_manifests = []
    summary = {
        "schema_version": 1,
        "representation": "boltz2 post-trunk s/z float16, recycle=1",
        "targets": {},
    }
    for target, config in TARGETS.items():
        target_dir = DESTINATION / target
        records_dir = target_dir / "records"
        provenance_dir = target_dir / "provenance"
        records_dir.mkdir(parents=True, exist_ok=True)
        provenance_dir.mkdir(parents=True, exist_ok=True)

        archive_provenance = []
        for archive_row in archives_by_target[target]:
            record = extract_archive(archive_row, records_dir)
            archive_provenance.append(record)
            archive_row["path"].unlink()
            archive_row["path"].with_suffix(".tar.sha256").unlink()

        for _index, source in sources[target].items():
            if source["kind"] == "run":
                move_local(source["path"], records_dir / source["path"].name)

        actual = sorted(records_dir.glob("embeddings_record_*.npz"))
        expected_names = {f"embeddings_record_{index:06d}.npz" for index in sources[target]}
        if {path.name for path in actual} != expected_names:
            raise ValueError(f"Final destination coverage mismatch: {target}")

        rows = []
        prep = preparation_manifest(target)
        prep_rows = {row["record_index"]: row for row in prep.iter_rows(named=True)}
        for index, source in sources[target].items():
            destination = records_dir / f"embeddings_record_{index:06d}.npz"
            row = prep_rows[index]
            digest = source.get("sha256") or sha256(destination)
            rows.append(
                {
                    "target": target,
                    "record_index": index,
                    "record_id": row["record_id"],
                    "compound_id": row["compound_id"],
                    "role": row["role"],
                    "standardized_smiles": row["standardized_smiles"],
                    "path": str(destination.relative_to(DESTINATION)),
                    "bytes": destination.stat().st_size,
                    "sha256": digest,
                    "source_kind": source["kind"],
                    "source_run_id": source.get("run_id"),
                    "source_archive": source["path"].name if source["kind"] == "archive" else None,
                }
            )
        manifest = pl.DataFrame(rows, infer_schema_length=None).sort("record_index")
        manifest.write_parquet(target_dir / "manifest.parquet")
        all_manifests.append(manifest)

        archive_dir = config["archives"]
        if archive_dir and archive_dir.exists():
            setup = archive_dir / "setup.json"
            if setup.exists():
                shutil.copy2(setup, provenance_dir / "colab-setup.json")
            for smoke in archive_dir.glob("smoke-*.tar*"):
                smoke.unlink()
            setup.unlink(missing_ok=True)
            archive_dir.rmdir()

        previous_summary_path = target_dir / "summary.json"
        previous_summary = (
            json.loads(previous_summary_path.read_text()) if previous_summary_path.exists() else {}
        )
        target_summary = {
            "records": manifest.height,
            "bytes": int(manifest.get_column("bytes").sum()),
            "preparation_run": config["preparation_run"],
            "source_runs": config["local_runs"],
            "source_archives": archive_provenance or previous_summary.get("source_archives", []),
            "manifest": str((target_dir / "manifest.parquet").relative_to(DESTINATION)),
        }
        (target_dir / "summary.json").write_text(json.dumps(target_summary, indent=2) + "\n")
        summary["targets"][target] = target_summary

    old_colab_root = ARTIFACTS / "boltz-trunk-colab"
    if old_colab_root.exists():
        old_colab_root.rmdir()
    combined = pl.concat(all_manifests, how="vertical_relaxed").sort("target", "record_index")
    combined.write_parquet(DESTINATION / "manifest.parquet")
    summary["records"] = combined.height
    summary["bytes"] = int(combined.get_column("bytes").sum())
    (DESTINATION / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (DESTINATION / "README.md").write_text(
        "# Boltz-2 post-trunk recycle-one tensors\n\n"
        "Canonical local store for CYP1A2, CYP2C9, CYP2D6 and CYP3A4. Each target's "
        "`records/` directory contains one `embeddings_record_XXXXXX.npz` per manifest "
        "row. `manifest.parquet` maps every file to its compound, role, source and hash.\n"
    )
    print(f"ORGANIZED records={combined.height} bytes={summary['bytes']} at {DESTINATION}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    organize(args.execute)


if __name__ == "__main__":
    main()
