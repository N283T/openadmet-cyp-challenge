"""Register the independent audit of selected generic-encoder D6 blends."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/regression"
inputs = [
    runs / "20260914T135408977818Z-a444a0b9/outputs/oof.parquet",
    runs / "20260912T083627802028Z-ae562a59/outputs/focus_2d6_combined_50pct.parquet",
    runs / "20260912T151747998390Z-f24622ac/outputs/crossfit-predictions.parquet",
]
for path in inputs:
    if not path.is_file():
        raise FileNotFoundError(path)

config = root / "configs/regression/d6-generic-encoder-screen-audit-v1.json"
payload = {
    "track": "regression",
    "endpoint": "CYP2D6_pIC50_direct_inhibition",
    "screen_run": "20260914T135408977818Z-a444a0b9",
    "expected_rows": 1457,
    "methods": [
        [view, method]
        for view in ["kermt", "molformer", "scores_kermt", "scores_molformer"]
        for method in ["pca64_ridge", "pca64_rbf_svr", "pls16", "extra_trees"]
    ],
    "candidates": [
        {
            "name": "molformer_extra_trees_match25",
            "view": "molformer",
            "method": "extra_trees",
            "placement": "match_latest_moments",
            "weight": 0.25,
        },
        {
            "name": "scores_molformer_extra_trees_match25",
            "view": "scores_molformer",
            "method": "extra_trees",
            "placement": "match_latest_moments",
            "weight": 0.25,
        },
        {
            "name": "scores_kermt_extra_trees_raw10",
            "view": "scores_kermt",
            "method": "extra_trees",
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
    "Audit generic KERMT and MoLFormer CYP2D6 blends by fold and molecular group",
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
