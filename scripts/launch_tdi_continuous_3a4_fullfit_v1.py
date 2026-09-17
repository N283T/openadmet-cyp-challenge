"""Launch test inference for both continuous-TDI fold-encoder components."""

import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "artifacts/runs"
SHARED = RUNS / "shared"
REG = RUNS / "regression"
common = {
    "measurements.parquet": SHARED / "20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    "assignments.parquet": REG
    / "20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
    "aliases.parquet": SHARED / "20260909T132734119233Z-a7063ccd/outputs/aliases.parquet",
    "reference.csv": ROOT / "data/raw/cyp-challenge-TEST-BLINDED.csv",
    "full-aux.parquet": SHARED / "20260912T074339967571Z-ff194fcf/outputs/full-aux.parquet",
}
families = {
    "arm": REG / "20260912T092823518869Z-8e57e8e7/outputs",
    "delta": REG / "20260912T092945291275Z-4fe962be/outputs",
}

for family, source in families.items():
    inputs = dict(common)
    inputs.update(
        {
            f"{family}-fold{fold}.parquet": source / f"fold{fold}/external-embeddings.parquet"
            for fold in range(5)
        }
    )
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        f"Full-label CYP3A4 readout over five {family} continuous-TDI fold encoders",
        "--config",
        str(
            ROOT / f"configs/submission-full/tdi_continuous_{family}_3a4_fold_encoder_lgbm_v1.json"
        ),
    ]
    for alias, path in inputs.items():
        command.extend(["--input", f"{alias}={path}"])
    result = subprocess.run(
        command + ["--", "python", "scripts/fullfit_fold_encoder_lgbm.py"],
        cwd=ROOT,
        text=True,
        capture_output=True,
    )
    print(result.stdout, flush=True)
    print(result.stderr, flush=True)
    result.check_returncode()
    run_ids = re.findall(r"run_id=(\S+)", result.stdout)
    assert len(run_ids) == 1
    print(f"{family}_run={run_ids[0]}", flush=True)
