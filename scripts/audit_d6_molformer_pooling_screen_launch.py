"""Register the independent audit of selected MoLFormer pooling blends."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/regression"
inputs = [
    runs / "20260914T140653173662Z-24d529d2/outputs/oof.parquet",
    runs / "20260912T083627802028Z-ae562a59/outputs/focus_2d6_combined_50pct.parquet",
    runs / "20260912T151747998390Z-f24622ac/outputs/crossfit-predictions.parquet",
]
for path in inputs:
    if not path.is_file():
        raise FileNotFoundError(path)

config = root / "configs/regression/d6-molformer-pooling-screen-audit-v1.json"
views = [
    "molformer_mean",
    "molformer_cls_mean",
    "scores_molformer_mean",
    "scores_molformer_cls_mean",
]
methods = [
    "extra_trees_l1",
    "extra_trees_l3",
    "extra_trees_l8",
    "extra_trees_f025",
    "random_forest",
    "catboost_mae",
    "xgboost",
]
payload = {
    "track": "regression",
    "endpoint": "CYP2D6_pIC50_direct_inhibition",
    "screen_run": "20260914T140653173662Z-24d529d2",
    "expected_rows": 1457,
    "methods": [[view, method] for view in views for method in methods],
    "candidates": [
        {
            "name": "scores_cls_mean_rf_match25",
            "view": "scores_molformer_cls_mean",
            "method": "random_forest",
            "placement": "match_latest_moments",
            "weight": 0.25,
        },
        {
            "name": "scores_cls_mean_rf_raw25",
            "view": "scores_molformer_cls_mean",
            "method": "random_forest",
            "placement": "raw",
            "weight": 0.25,
        },
        {
            "name": "scores_cls_mean_extra_trees_l1_raw25",
            "view": "scores_molformer_cls_mean",
            "method": "extra_trees_l1",
            "placement": "raw",
            "weight": 0.25,
        },
        {
            "name": "scores_cls_mean_catboost_match25",
            "view": "scores_molformer_cls_mean",
            "method": "catboost_mae",
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
    "Audit selected MoLFormer pooling CYP2D6 blends by fold and molecular group",
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
