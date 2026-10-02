"""Declare every own model weight and input for the concrete candidate review."""

import subprocess
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    runs = root / "artifacts/runs"
    submissions = root / "artifacts/submissions/classification"
    output = root / "artifacts/runs/classification/20261002T090029588080Z-2a20da72/outputs"
    inputs = {
        "blind-policy.parquet": output / "blind-policy.parquet",
        "blind-folds.parquet": output / "blind-fold-predictions.parquet",
        "features.parquet": runs
        / "shared/20260912T061231617816Z-299439e8/outputs/external-embeddings.parquet",
        "last-sent.parquet": submissions / "20261001T141041310215Z-c7faa79f/predictions.parquet",
        "best-control.parquet": submissions / "20260924T050705396203Z-10e1d3ad/predictions.parquet",
        "reference.csv": root / "data/raw/cyp-challenge-TEST-BLINDED.csv",
        "encoder-provenance.json": root
        / "artifacts/runs/shared/20260912T061231617816Z-299439e8/config.json",
    }
    for seed in [20261002, 20261017]:
        for bag in range(5):
            inputs[f"model-{seed}-bag{bag}.pt"] = output / f"models/{seed}-bag{bag}.pt"
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        "Verify ten persisted new readouts and freeze a concrete D6 candidate with exact "
        "A4 retention and target-wise change counts; no upload",
        "--config",
        str(root / "configs/classification/d6-chemeleon-candidate-package-v1.json"),
    ]
    for alias, path in inputs.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        command.extend(["--input", f"{alias}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/package_cyp_chemeleon_candidate.py"],
        cwd=root,
        check=True,
    )


if __name__ == "__main__":
    main()
