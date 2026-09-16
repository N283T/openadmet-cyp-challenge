"""Register the independent audit of selected CYP2D6 pairwise-ranker blends."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/regression"
inputs = [
    runs / "20260914T123940860032Z-de101db4/outputs/oof.parquet",
    runs / "20260912T083627802028Z-ae562a59/outputs/focus_2d6_combined_50pct.parquet",
    runs / "20260912T151747998390Z-f24622ac/outputs/crossfit-predictions.parquet",
]
for path in inputs:
    if not path.is_file():
        raise FileNotFoundError(path)

config = root / "configs/regression/d6-pairwise-ranker-screen-audit-v1.json"
payload = {
    "track": "regression",
    "endpoint": "CYP2D6_pIC50_direct_inhibition",
    "screen_run": "20260914T123940860032Z-de101db4",
    "expected_rows": 1457,
    "methods": [
        [view, method]
        for view in ["scores51", "medium604", "hybrid2699"]
        for method in [
            "pairwise_logistic",
            "rank_svm",
            "ranknet",
            "xgb_pairwise",
            "cat_pairlogit",
        ]
    ],
    "candidates": [
        {
            "name": "scores_xgb_pairwise_10",
            "view": "scores51",
            "method": "xgb_pairwise",
            "placement": "raw",
            "weight": 0.1,
        },
        {
            "name": "scores_pairwise_logistic_10",
            "view": "scores51",
            "method": "pairwise_logistic",
            "placement": "raw",
            "weight": 0.1,
        },
        {
            "name": "scores_pairwise_logistic_25",
            "view": "scores51",
            "method": "pairwise_logistic",
            "placement": "raw",
            "weight": 0.25,
        },
        {
            "name": "hybrid_cat_pairlogit_10",
            "view": "hybrid2699",
            "method": "cat_pairlogit",
            "placement": "raw",
            "weight": 0.1,
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
    "Audit selected CYP2D6 pairwise-ranker blends across folds and molecular groups",
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
