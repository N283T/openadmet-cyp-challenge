"""Record source-interval and regime diagnostics for the independent restart."""

import subprocess
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    pilot = root / "artifacts/runs/classification/20261002T004016766672Z-e2d420d4"
    inputs = {
        name: pilot / "outputs" / name
        for name in [
            "summary.json",
            "metrics.csv",
            "operating-points.csv",
            "oof.parquet",
            "assignments.parquet",
            "source-audit.parquet",
            "representation.json",
        ]
    }
    inputs["measurements.parquet"] = pilot / "inputs/measurements.parquet"
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        "Diagnose independent D6 model uncertainty and label-regime failures; no old model anchor",
        "--config",
        str(root / "configs/classification/d6-zero-start-diagnosis-v1.json"),
    ]
    for alias, path in inputs.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        command.extend(["--input", f"{alias}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/diagnose_cyp_d6_zero_start.py"], cwd=root, check=True
    )


if __name__ == "__main__":
    main()
