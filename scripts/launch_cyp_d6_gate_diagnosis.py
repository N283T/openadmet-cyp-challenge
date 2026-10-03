"""Launch the read-only diagnosis of the D6 activity gate on blind chemistry."""

import argparse
import subprocess
from pathlib import Path

DATA = "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
SPLITS = "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs"
POOL_RUN = "artifacts/runs/classification/20261002T215345688167Z-7920b49a/outputs"
TESTLIKE = "artifacts/tracked/cyp-d6-testlike-pool-20261002T215345688167Z-7920b49a"
INPUTS = {
    "testlike-pool-oof.parquet": f"{TESTLIKE}/oof_scores.parquet",
    "test_scores.parquet": f"{POOL_RUN}/test_scores.parquet",
    "aliases.parquet": f"{DATA}/aliases.parquet",
    "primary-screen.csv": "data/raw/cyp-challenge-single-concentration-TRAIN.csv",
    "assignments.parquet": f"{SPLITS}/assignments.parquet",
}
PURPOSE = (
    "Diagnose the D6 activity gate before trusting or discarding it as a call filter: "
    "blind versus development p_active distributions by primary-screen stratum, held-fold "
    "reliability of the gate against observed primary activity, the unscreened development "
    "stratum, the officially labelled positive rate by p_active bin, a primary-region "
    "decomposition giving an implied blind positive rate, and the circularity of a "
    "pseudo-negative-labelled metric for judging a gate filter; read-only, no candidate or upload"
)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/classification/d6-gate-diagnosis-v1.json")
    args = parser.parse_args()
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        PURPOSE,
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
        command + ["--", "python", "scripts/evaluate_cyp_d6_gate_diagnosis.py"],
        cwd=root,
        check=True,
    )


if __name__ == "__main__":
    main()
