"""Register and launch combinations of retained CYP2D6 ordering signals."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/regression"
inputs = [
    runs / "20260912T083627802028Z-ae562a59/outputs/focus_2d6_combined_50pct.parquet",
    runs / "20260912T151747998390Z-f24622ac/outputs/crossfit-predictions.parquet",
    runs / "20260914T140653173662Z-24d529d2/outputs/oof.parquet",
    runs / "20260914T121807545564Z-5d24aefa/outputs/oof.parquet",
    runs / "20260914T125646427763Z-1707a0e2/outputs/oof.parquet",
    runs / "20260914T123940860032Z-de101db4/outputs/oof.parquet",
]
for path in inputs:
    if not path.is_file():
        raise FileNotFoundError(path)

config = root / "configs/regression/d6-candidate-stack-v1.json"
payload = {
    "track": "regression",
    "endpoint": "CYP2D6_pIC50_direct_inhibition",
    "expected_rows": 1457,
    "source_runs": [
        "20260914T140653173662Z-24d529d2",
        "20260914T121807545564Z-5d24aefa",
        "20260914T125646427763Z-1707a0e2",
        "20260914T123940860032Z-de101db4",
    ],
    "fixed_components": [
        "MoLFormer RF25 moment matched",
        "hybrid CatBoost50 moment matched",
        "ordinal-logit5 raw10",
        "scores XGBoost-pairwise raw10",
    ],
    "meta_models": ["Ridge alpha 1/10/100", "positive linear regression"],
    "bootstrap_replicates": 5000,
    "bootstrap_seed": 20260914,
    "evaluation": "Fixed combinations plus five-fold-excluded level-2 stacking.",
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
    "Combine retained MoLFormer, CatBoost, ordinal, and pairwise CYP2D6 signals",
    "--config",
    str(config),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/regression_d6_candidate_stack.py"],
    cwd=root,
    check=True,
)
