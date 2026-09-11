"""Launch the predeclared five-fold external-data comparison through durable runs."""

import argparse
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--folds", nargs="+", type=int, default=list(range(5)))
    parser.add_argument(
        "--variants", nargs="+", choices=["control", "external"], default=["control", "external"]
    )
    args = parser.parse_args()
    base = json.loads((ROOT / "configs/external/training-v1.json").read_text())
    runs = ROOT / "artifacts/runs/shared"
    inputs = [
        runs / base["dataset_run"] / "outputs" / name
        for name in ["compounds.parquet", "measurements.parquet"]
    ]
    inputs += [
        runs / base["split_run"] / "outputs/assignments.parquet",
        runs / base["auxiliary_run"] / "outputs/targets.parquet",
    ]
    inputs += [
        runs / base["external_run"] / "outputs" / name
        for name in [
            "external-compounds.parquet",
            "external-targets.parquet",
            "challenge-connectivity.parquet",
        ]
    ]
    inputs += [runs / base["checkpoint_run"] / "outputs/chemeleon_mp.pt"]
    for fold in args.folds:
        for variant in args.variants:
            cfg = {
                **base,
                "fold": fold,
                "variant": variant,
                "external_weight": 0.3 if variant == "external" else 0.0,
            }
            path = ROOT / f"configs/external/{variant}-fold{fold}-v1.json"
            if path.exists():
                assert json.loads(path.read_text()) == cfg
            else:
                path.write_text(json.dumps(cfg, indent=2) + "\n")
            command = [
                "pixi",
                "run",
                "cyp",
                "run",
                "--track",
                "regression",
                "--purpose",
                f"External-data comparison {variant}, UMAP fold {fold}; ChEMBL/PubChem auxiliaries",
                "--config",
                str(path),
            ]
            for file in inputs:
                command += ["--input", str(file)]
            subprocess.run(
                command + ["--", "python", "scripts/external_regression.py"], cwd=ROOT, check=True
            )


if __name__ == "__main__":
    main()
