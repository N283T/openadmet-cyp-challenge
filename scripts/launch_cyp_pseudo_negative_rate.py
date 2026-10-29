"""Launch the direct, non-circular test of the pseudo-negative assumption."""

import argparse
import subprocess
from pathlib import Path

DATA = "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
TESTLIKE = "artifacts/tracked/cyp-d6-testlike-pool-20261002T215345688167Z-7920b49a"
INPUTS = {
    "compounds.parquet": f"{DATA}/compounds.parquet",
    "measurements.parquet": f"{DATA}/measurements.parquet",
    "aliases.parquet": f"{DATA}/aliases.parquet",
    "test_mapping.parquet": f"{DATA}/test_mapping.parquet",
    "primary-screen.csv": "data/raw/cyp-challenge-single-concentration-TRAIN.csv",
    "emax-train.csv": "data/raw/cyp-challenge-TRAIN_Emax.csv",
    "testlike-pool-oof.parquet": f"{TESTLIKE}/oof_scores.parquet",
}


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/classification/pseudo-negative-rate-v1.json")
    args = parser.parse_args()
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        "Pseudo-negative assumption, measured directly: positive rate among labelled "
        "weak-primary 1A2/2C9/3A4 rows with exact binomial bounds and a log2fc cut sweep, "
        "plus the number of hidden D6 positives needed to overturn the pool-versus-official "
        "ranking; nothing fitted, no upload",
        "--config",
        str(root / args.config),
    ]
    for alias, relative in INPUTS.items():
        path = root / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        command.extend(["--input", f"{alias}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/evaluate_cyp_pseudo_negative_rate.py"],
        cwd=root,
        check=True,
    )


if __name__ == "__main__":
    main()
