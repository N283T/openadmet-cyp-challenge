"""Record missing-outcome identification bounds and their controlled A4 check."""

import subprocess
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    runs = root / "artifacts/runs/classification"
    inputs = {}
    for alias, run, filename in [
        ("diverse-oof.parquet", "20260913T112957931554Z-32d00f6e", "oof.parquet"),
        ("fixed-oof.parquet", "20260913T123045515718Z-0df22d54", "oof.parquet"),
        ("control-umap.parquet", "20260924T040301085704Z-914bb06c", "all-pool-oof.parquet"),
        ("control-scaffold.parquet", "20260924T041452692950Z-e17778e3", "all-pool-oof.parquet"),
        ("ap-umap.parquet", "20260924T035023873825Z-b719ba35", "all-pool-oof.parquet"),
        ("ap-scaffold.parquet", "20260924T041348966036Z-8db19b4b", "all-pool-oof.parquet"),
        ("a4-oof.parquet", "20261002T042725116238Z-fe96ab3c", "nested-oof.parquet"),
        ("a4-mask.parquet", "20261002T042725116238Z-fe96ab3c", "selection-control-source.parquet"),
    ]:
        inputs[alias] = runs / run / "outputs" / filename
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        "Calculate exact feasible MCC-difference bounds under missing classes; "
        "test A4 population-reversal coverage and historical full-pool comparability",
        "--config",
        str(root / "configs/classification/d6-partial-label-bounds-v1.json"),
    ]
    for alias, path in inputs.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        command.extend(["--input", f"{alias}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/evaluate_cyp_partial_label_bounds.py"],
        cwd=root,
        check=True,
    )


if __name__ == "__main__":
    main()
