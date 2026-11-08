"""Copy saved leaderboard CSVs into a tracked snapshot directory.

Each argument after the output directory is a completed `cyp leaderboard` API
directory under `artifacts/api/`. The CSV is copied byte for byte to
`csv/<board>/<track>/<endpoint>.csv` and checked against the SHA-256 recorded in
the API manifest. `manifest.json` and `tracked-files.json` follow the format of
`artifacts/tracked/leaderboard-snapshot-20260926/`.
"""

import argparse
import hashlib
import json
import shutil
from datetime import UTC, datetime
from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("out", type=Path)
    parser.add_argument("--captured-on", required=True)
    parser.add_argument("--scope", required=True)
    parser.add_argument("api_dirs", nargs="+", type=Path)
    args = parser.parse_args()

    out = args.out if args.out.is_absolute() else ROOT / args.out
    if out.exists():
        raise SystemExit(f"{out} already exists; snapshots are immutable")
    entries, files = [], []
    for api_dir in args.api_dirs:
        api_dir = (ROOT / api_dir).resolve() if not api_dir.is_absolute() else api_dir
        manifest = json.loads((api_dir / "manifest.json").read_text())
        if manifest.get("status") != "completed":
            raise SystemExit(f"{api_dir} is not a completed download")
        csvs = [f for f in manifest["files"] if f["path"].endswith(".csv")]
        if len(csvs) != 1:
            raise SystemExit(f"{api_dir} has {len(csvs)} CSV files")
        record = csvs[0]
        source = api_dir / record["path"]
        if sha256(source) != record["sha256"]:
            raise SystemExit(f"{source} does not match its manifest SHA-256")
        relative_source = str(source.relative_to(ROOT))
        target = out / "csv" / manifest["board"] / manifest["track"] / f"{manifest['endpoint']}.csv"
        if target.exists():
            raise SystemExit(f"duplicate board/track/endpoint: {target}")
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        if sha256(target) != record["sha256"] or target.stat().st_size != record["bytes"]:
            raise SystemExit(f"{target} copy check failed")
        entries.append(
            {
                "board": manifest["board"],
                "track": manifest["track"],
                "endpoint": manifest["endpoint"],
                "captured_at": manifest["created_at"],
                "api_name": manifest["api_name"],
                "csv": relative_source,
                "sha256": record["sha256"],
                "bytes": record["bytes"],
                "rows": pl.read_csv(source, infer_schema_length=0).height,
            }
        )
        files.append(
            {
                "board": manifest["board"],
                "track": manifest["track"],
                "endpoint": manifest["endpoint"],
                "path": str(target.relative_to(out)),
                "sha256": record["sha256"],
                "source_api_path": relative_source,
            }
        )
    key = lambda item: (item["board"], item["track"], item["endpoint"])  # noqa: E731
    manifest = {
        "captured_on": args.captured_on,
        "indexed_at": datetime.now(UTC).isoformat(),
        "source": "https://huggingface.co/spaces/openadmet/cyp-challenge",
        "scope": args.scope,
        "note": (
            "Raw public leaderboard CSVs and API response/config snapshots are kept locally "
            "under artifacts/api/ and are Git-ignored. Paths are relative to the repository root."
        ),
        "entries": sorted(entries, key=key),
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    tracked = {"source_manifest": "manifest.json", "files": sorted(files, key=key)}
    (out / "tracked-files.json").write_text(json.dumps(tracked, indent=2) + "\n")
    print(out)


if __name__ == "__main__":
    main()
