"""Promote the pseudo-negative validation outputs without overwriting an existing record."""

import argparse
import hashlib
import json
import shutil
from pathlib import Path

PREFIXES = {
    "evaluate_cyp_pseudo_negative_rate.py": "cyp-pseudo-negative-rate",
    "evaluate_cyp_induced_defect.py": "cyp-induced-defect",
    "evaluate_cyp_d6_official_metric.py": "cyp-d6-official-metric",
    "evaluate_cyp_d6_probe_preregistration.py": "cyp-d6-probe-preregistration",
}


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
        manifest = json.loads((run / "manifest.json").read_text())
        script = next((Path(part).name for part in manifest["command"] if part.endswith(".py")), "")
        prefix = PREFIXES.get(script)
        if prefix is None:
            raise KeyError(f"Unknown launcher for retention: {script!r} ({run_id})")
        target = root / "artifacts/tracked" / f"{prefix}-{run_id}"
        if target.exists():
            raise FileExistsError(target)
        shutil.copytree(run / "outputs", target)
        hashes = {
            path.relative_to(target).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(target.rglob("*"))
            if path.is_file()
        }
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
