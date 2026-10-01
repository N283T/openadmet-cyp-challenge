"""Record multi-trunk placement and pair correlations on incumbent cohorts."""

import hashlib
import json
import shutil
from pathlib import Path

from launch_cyp_boltz_incumbent_suite import recorded


def main():
    root = Path(__file__).resolve().parents[1]
    suite = root / "artifacts/boltz-trunk/experiments/20261001T065518Z-incumbent"
    source_id = (suite / "regression-placement-run.txt").read_text().strip()
    source = root / "artifacts/runs/regression" / source_id
    cfg = json.loads((source / "config.json").read_text())
    cfg["block_sizes"] = [1, 2, 3, 5, 10]
    cfg["blend_weights"] = [0, 0.05, 0.1, 0.2, 0.3, 0.5, 1]
    cfg["source_placement_run"] = source_id
    inputs = {
        Path(item["snapshot"]).name: Path(item["snapshot"])
        for item in json.loads((source / "manifest.json").read_text())["inputs"]
    }
    identifier = recorded(
        suite, "multi-trunk-ensemble", "regression", cfg, inputs, "evaluate_cyp_boltz_multi.py"
    )
    (suite / "multi-trunk-run.txt").write_text(identifier + "\n")
    output = root / "artifacts/runs/regression" / identifier / "outputs"
    tracked = root / "artifacts/tracked" / f"cyp-boltz-multi-{identifier}"
    tracked.mkdir()
    inventory = []
    for path in output.iterdir():
        shutil.copy2(path, tracked / path.name)
        inventory.append(
            {"file": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        )
    (tracked / "provenance.json").write_text(
        json.dumps(
            {"run": identifier, "source_placement_run": source_id, "files": inventory}, indent=2
        )
        + "\n"
    )
    print(f"multi_trunk_run={identifier}")


if __name__ == "__main__":
    main()
