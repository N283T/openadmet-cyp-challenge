"""Launch an independent molecular/measurement D6 pilot from canonical source data."""

import argparse
import json
import subprocess
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/classification/d6-zero-start-v1.json")
    args = parser.parse_args()
    data = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        "Restart D6 TDI from source molecules and paired measurements: independent fresh "
        "models, nested thresholds with no call budget, random and structural holdouts",
        "--config",
        str(root / args.config),
    ]
    for name in ["compounds.parquet", "measurements.parquet", "aliases.parquet"]:
        path = data / name
        if not path.is_file():
            raise FileNotFoundError(path)
        command.extend(["--input", f"{name}={path}"])
    cfg = json.loads((root / args.config).read_text())
    if cfg.get("fixed_assignments"):
        path = (
            root
            / "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet"
        )
        command.extend(["--input", f"{cfg['fixed_assignments']}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/evaluate_cyp_d6_zero_start.py"], cwd=root, check=True
    )


if __name__ == "__main__":
    main()
