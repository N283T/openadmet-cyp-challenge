"""Record the full candidate decision algorithm on untouched outer groups."""

import argparse
import json
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config", default="configs/classification/d6-nested-bag-validation-v1.json"
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    runs = root / "artifacts/runs"
    shared = runs / "shared/20260909T132734119233Z-a7063ccd/outputs"
    inputs = {
        name: shared / name
        for name in ["compounds.parquet", "measurements.parquet", "aliases.parquet"]
    }
    inputs["assignments.parquet"] = (
        runs / "shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet"
    )
    cfg = json.loads((root / args.config).read_text())
    if cfg.get("primary_screen"):
        inputs[cfg["primary_screen"]] = (
            root / "data/raw/cyp-challenge-single-concentration-TRAIN.csv"
        )
    if cfg.get("embedding_file"):
        inputs[cfg["embedding_file"]] = (
            runs / "shared/20260912T061231617816Z-299439e8/outputs/external-embeddings.parquet"
        )
    for alias, run in {
        "control-umap.parquet": "20260924T040301085704Z-914bb06c",
        "control-scaffold.parquet": "20260924T041452692950Z-e17778e3",
        "ap-umap.parquet": "20260924T035023873825Z-b719ba35",
        "ap-scaffold.parquet": "20260924T041348966036Z-8db19b4b",
    }.items():
        if alias in cfg.get("historical_comparisons", []):
            inputs[alias] = runs / "classification" / run / "outputs/all-pool-oof.parquet"
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        cfg.get("scientific_purpose", "Validate the full nested decision pipeline"),
        "--config",
        str(root / args.config),
    ]
    for alias, path in inputs.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        command.extend(["--input", f"{alias}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/evaluate_cyp_d6_zero_start.py"], cwd=root, check=True
    )


if __name__ == "__main__":
    main()
