"""Promote the unsubmitted blind review and exact local candidate identities."""

import argparse
import hashlib
import json
import shutil
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("package_run")
    parser.add_argument("candidate_id")
    parser.add_argument("prepared_id")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    runs = root / "artifacts/runs/classification"
    source = runs / args.package_run
    if json.loads((source / "result.json").read_text())["status"] != "completed":
        raise ValueError("Package run is incomplete")
    summary = json.loads((source / "outputs/summary.json").read_text())
    selection = runs / summary["selection_run"]
    dry = json.loads((source / "scratch/prepared-dry-run.json").read_text())
    if dry["sent"] or dry["submission_id"] != args.prepared_id:
        raise ValueError("Expected the exact unsent dry-run receipt")
    manifests = {}
    for kind, identifier in [("candidates", args.candidate_id), ("submissions", args.prepared_id)]:
        path = root / "artifacts" / kind / "classification" / identifier / "manifest.json"
        manifest = json.loads(path.read_text())
        if manifest["sha256"] != summary["prediction_sha256"]:
            raise ValueError("Frozen payload hash differs")
        manifests[kind] = manifest
    tracked = root / "artifacts/tracked" / f"cyp-boltz-d6-ap-probe-{args.package_run}"
    tracked.mkdir()
    for path in (source / "outputs").iterdir():
        if path.is_file():
            shutil.copy2(path, tracked / path.name)
    for path in (selection / "outputs").iterdir():
        if path.is_file():
            shutil.copy2(path, tracked / ("selection-" + path.name))
    (tracked / "prepared-dry-run.json").write_text(json.dumps(dry, indent=2) + "\n")
    provenance = {
        "package_run": args.package_run,
        "selection_run": summary["selection_run"],
        "candidate_id": args.candidate_id,
        "prepared_id": args.prepared_id,
        "uploaded": False,
        "frozen_manifests": manifests,
        "package_inputs": json.loads((source / "manifest.json").read_text())["inputs"],
        "files": {
            path.name: {
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "bytes": path.stat().st_size,
            }
            for path in sorted(tracked.iterdir())
            if path.is_file()
        },
    }
    (tracked / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(tracked)


if __name__ == "__main__":
    main()
