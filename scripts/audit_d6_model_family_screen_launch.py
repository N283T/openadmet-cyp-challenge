"""Register the independent audit of selected CYP2D6 model-family screen rows."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/regression"
inputs = [
    runs / "20260914T121807545564Z-5d24aefa/outputs/oof.parquet",
    runs / "20260912T083627802028Z-ae562a59/outputs/focus_2d6_combined_50pct.parquet",
    runs / "20260912T151747998390Z-f24622ac/outputs/crossfit-predictions.parquet",
]
for path in inputs:
    if not path.is_file():
        raise FileNotFoundError(path)

config = root / "configs/regression/d6-model-family-screen-audit-v1.json"
payload = {
    "track": "regression",
    "endpoint": "CYP2D6_pIC50_direct_inhibition",
    "screen_run": "20260914T121807545564Z-5d24aefa",
    "expected_rows": 1457,
    "methods": [
        [view, method]
        for view in ["scores51", "medium604", "hybrid2699"]
        for method in [
            "ridge",
            "elastic_net",
            "pls16",
            "extra_trees",
            "random_forest",
            "rbf_svr",
            "nystroem_ridge",
            "knn",
            "mlp",
            "xgboost",
            "catboost",
            "compact_kan",
        ]
    ],
    "candidates": [
        {
            "name": "hybrid_catboost_50",
            "view": "hybrid2699",
            "method": "catboost",
            "placement": "match_latest_moments",
            "weight": 0.5,
        },
        {
            "name": "scores_random_forest_25",
            "view": "scores51",
            "method": "random_forest",
            "placement": "match_latest_moments",
            "weight": 0.25,
        },
        {
            "name": "scores_pls16_25",
            "view": "scores51",
            "method": "pls16",
            "placement": "match_latest_moments",
            "weight": 0.25,
        },
        {
            "name": "medium_catboost_25",
            "view": "medium604",
            "method": "catboost",
            "placement": "match_latest_moments",
            "weight": 0.25,
        },
    ],
    "bootstrap_replicates": 5000,
    "bootstrap_seed": 20260914,
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
    "Audit selected CYP2D6 model-family blends across folds and molecular-group bootstrap",
    "--config",
    str(config),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/audit_d6_model_family_screen.py"],
    cwd=root,
    check=True,
)
