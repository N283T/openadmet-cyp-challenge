"""Launch the D6 test-like pool evaluation (primary-screen pseudo-negatives + gate)."""

import argparse
import subprocess
from pathlib import Path

DATA = "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
SPLITS = "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs"
CHEMELEON = "artifacts/runs/shared/20260910T012630927160Z-e9576045/outputs"
SEP15 = "artifacts/submissions/classification/20260924T050705396203Z-10e1d3ad"
LEGACY = "artifacts/runs/classification/20261002T120413515280Z-f2330789/outputs"
INPUTS = {
    "compounds.parquet": f"{DATA}/compounds.parquet",
    "measurements.parquet": f"{DATA}/measurements.parquet",
    "aliases.parquet": f"{DATA}/aliases.parquet",
    "test_mapping.parquet": f"{DATA}/test_mapping.parquet",
    "assignments.parquet": f"{SPLITS}/assignments.parquet",
    "primary-screen.csv": "data/raw/cyp-challenge-single-concentration-TRAIN.csv",
    "chemeleon-embeddings.parquet": f"{CHEMELEON}/embeddings.parquet",
    "sep15-d6-payload.parquet": f"{SEP15}/predictions.parquet",
    "legacy-sep15-nested-oof.parquet": f"{LEGACY}/nested-oof.parquet",
}


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/classification/d6-testlike-pool-v1.json")
    args = parser.parse_args()
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        "D6 TDI test-like pool: add primary-screen non-hits (log2fc >= -0.25) as "
        "pseudo-negatives, compare official-only, pool-trained and activity-gated "
        "models on umap/scaffold complete pools with nested rate selection; no upload",
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
        command + ["--", "python", "scripts/evaluate_cyp_d6_testlike_pool.py"],
        cwd=root,
        check=True,
    )


if __name__ == "__main__":
    main()
