"""Record the completed, matched primary-selection and potency-gate analyses."""

import json
import subprocess
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    config = root / "configs/classification/primary-selection-controls-review-v1.json"
    cfg = json.loads(config.read_text())
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        "Measure selected-population versus broad-population reversal on known A4 "
        "outcomes and compare D6 primary and assay-compatible potency gates",
        "--config",
        str(config),
    ]
    for case in cfg["cases"]:
        out = root / "artifacts/runs/classification" / case["run_id"] / "outputs"
        for key, filename in [
            ("predictions", "nested-oof.parquet"),
            (
                "source",
                "selection-control-source.parquet"
                if case["selection_control"]
                else "source-audit.parquet",
            ),
        ]:
            path = out / filename
            if not path.is_file():
                raise FileNotFoundError(path)
            command.extend(["--input", f"{case[key]}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/report_cyp_primary_selection_controls.py"],
        cwd=root,
        check=True,
    )


if __name__ == "__main__":
    main()
