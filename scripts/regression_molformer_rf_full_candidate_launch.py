"""Register the full-fit MoLFormer RF D6 candidate run."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = [
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    root / "artifacts/runs/shared/20260912T072105522912Z-53eb675c/outputs/external-scores.parquet",
    root / "artifacts/runs/regression/20260914T140517644424Z-de551457/outputs/features.parquet",
    root / "artifacts/runs/shared/20260912T074339967571Z-ff194fcf/outputs/full-aux.parquet",
    root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/aliases.parquet",
    root / "artifacts/candidates/regression/20260914T024735834387Z-f87e9914/predictions.csv",
    root / "data/raw/cyp-challenge-TEST-BLINDED.csv",
]
for path in inputs:
    if not path.is_file():
        raise FileNotFoundError(path)

config = root / "configs/regression/d6-molformer-rf-full-candidate-v1.json"
payload = {
    "track": "regression",
    "endpoint": "CYP2D6_pIC50_direct_inhibition",
    "selection_run": "20260914T140653173662Z-24d529d2",
    "audit_run": "20260914T142335850035Z-781cf3a4",
    "encoder_run": "20260914T140517644424Z-de551457",
    "base_candidate": "20260914T024735834387Z-f87e9914",
    "expected_train_rows": 1493,
    "tree_estimators": 240,
    "min_samples_leaf": 3,
    "max_features": 0.5,
    "threads": 4,
    "seed": 20260914,
    "blend_weight": 0.25,
    "placement": "Match full-fit RF test moments to the live D6 vector before blending.",
    "other_endpoints": "Bitwise unchanged from the currently scored regression candidate.",
}
config.write_text(json.dumps(payload, indent=2) + "\n")

command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Full-fit and package the selected MoLFormer Random Forest D6 25% blend",
    "--config",
    str(config),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/regression_molformer_rf_full_candidate.py"],
    cwd=root,
    check=True,
)
