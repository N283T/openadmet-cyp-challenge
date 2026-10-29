"""Launch the induced D6 coverage defect and pseudo-negative fix experiment."""

import argparse
import subprocess
from pathlib import Path

DATA = "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
SPLITS = "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs"
INPUTS = {
    "compounds.parquet": f"{DATA}/compounds.parquet",
    "measurements.parquet": f"{DATA}/measurements.parquet",
    "aliases.parquet": f"{DATA}/aliases.parquet",
    "test_mapping.parquet": f"{DATA}/test_mapping.parquet",
    "assignments.parquet": f"{SPLITS}/assignments.parquet",
    "primary-screen.csv": "data/raw/cyp-challenge-single-concentration-TRAIN.csv",
}


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/classification/induced-defect-v1.json")
    args = parser.parse_args()
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        "Induce the D6 label-coverage defect on 1A2/2C9/3A4 by hiding weak and "
        "missing-primary labels from fitting and rate selection, then re-add the "
        "weak-primary rows as pseudo-negatives and judge both arms plus a "
        "real-complete-label ceiling on the real held-fold labels; no upload",
        "--config",
        str(root / args.config),
        "--input-mode",
        "hardlink",
    ]
    for alias, relative in INPUTS.items():
        path = root / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        command.extend(["--input", f"{alias}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/evaluate_cyp_induced_defect.py"],
        cwd=root,
        check=True,
    )


if __name__ == "__main__":
    main()
