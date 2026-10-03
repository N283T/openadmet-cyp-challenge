"""Launch the D6 comparison under metrics that do not score pseudo-negative labels."""

import argparse
import subprocess
from pathlib import Path

SPLITS = "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs"
TESTLIKE = "artifacts/tracked/cyp-d6-testlike-pool-20261002T215345688167Z-7920b49a"
INPUTS = {
    "testlike-pool-oof.parquet": f"{TESTLIKE}/oof_scores.parquet",
    "assignments.parquet": f"{SPLITS}/assignments.parquet",
}


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/classification/d6-official-metric-v1.json")
    args = parser.parse_args()
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        "Rescore the frozen D6 held-fold predictions under an official-labels-only "
        "prevalence-reweighted metric and a stratified mixture metric whose mixture "
        "fraction is identified from the published prevalence, on a rate grid extended "
        "to 50%; nothing fitted, no upload",
        "--config",
        str(root / args.config),
    ]
    for alias, relative in INPUTS.items():
        path = root / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        command.extend(["--input", f"{alias}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/evaluate_cyp_d6_official_metric.py"],
        cwd=root,
        check=True,
    )


if __name__ == "__main__":
    main()
