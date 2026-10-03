"""Promote the D6 Live gain post-mortem outputs without overwriting an existing record."""

import argparse
import hashlib
import json
import shutil
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_ids", nargs="+")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    for run_id in args.run_ids:
        run = root / "artifacts/runs/classification" / run_id
        result = json.loads((run / "result.json").read_text())
        if result.get("exit_code") != 0:
            raise ValueError(f"Cannot retain incomplete run: {run_id}")
        target = root / "artifacts/tracked" / f"cyp-d6-live-gain-{run_id}"
        if target.exists():
            raise FileExistsError(target)
        shutil.copytree(run / "outputs", target)
        hashes = {
            path.relative_to(target).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(target.rglob("*"))
            if path.is_file()
        }
        manifest = json.loads((run / "manifest.json").read_text())
        provenance = dict(
            source_run=f"classification/{run_id}",
            purpose=manifest.get("purpose") or manifest.get("intent", {}).get("purpose"),
            outputs=hashes,
            input_provenance=manifest["inputs"],
            config=json.loads((run / "config.json").read_text()),
        )
        (target / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
        print(target)


if __name__ == "__main__":
    main()
