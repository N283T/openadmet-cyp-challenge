"""Record selected trunk/incumbent correlation without any model refits."""

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
    cfg["targets"] = {t: cfg["targets"][t] for t in ["cyp1a2", "cyp2d6"]}
    inputs = {}
    for item in json.loads((source / "manifest.json").read_text())["inputs"]:
        alias = Path(item["snapshot"]).name
        if alias.startswith(("cyp1a2-", "cyp2d6-", "d6-base")):
            inputs[alias] = Path(item["snapshot"])
    inputs["incumbent-oof.parquet"] = source / "outputs/oof.parquet"
    inputs["selections.parquet"] = source / "outputs/selections.parquet"
    identifier = recorded(
        suite,
        "selected-trunk-incumbent-correlations",
        "regression",
        cfg,
        inputs,
        "cyp_boltz_incumbent_correlation.py",
    )
    (suite / "correlation-run.txt").write_text(identifier + "\n")
    output = root / "artifacts/runs/regression" / identifier / "outputs"
    tracked = root / "artifacts/tracked" / f"cyp-boltz-correlations-{identifier}"
    tracked.mkdir()
    inventory = []
    for path in output.iterdir():
        shutil.copy2(path, tracked / path.name)
        inventory.append(
            {"file": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        )
    (tracked / "provenance.json").write_text(
        json.dumps(
            {
                "run": identifier,
                "source_placement_run": source_id,
                "files": inventory,
            },
            indent=2,
        )
        + "\n"
    )
    print(f"correlation_run={identifier}")


if __name__ == "__main__":
    main()
