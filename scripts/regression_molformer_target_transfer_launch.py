"""Apply the D6-selected MoLFormer view to the other regression endpoints."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
inputs = [
    runs / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
    runs / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    runs / "shared/20260912T072105522912Z-53eb675c/outputs/external-scores.parquet",
    runs / "shared/20260912T082554377951Z-d272774c/outputs/nested-log2fc.parquet",
    runs / "regression/20260912T083627802028Z-ae562a59/outputs/focus_2d6_combined_50pct.parquet",
    runs / "regression/20260912T151747998390Z-f24622ac/outputs/crossfit-predictions.parquet",
    runs / "regression/20260914T140517644424Z-de551457/outputs/features.parquet",
]
for path in inputs:
    if not path.is_file():
        raise FileNotFoundError(path)

endpoints = [
    "CYP1A2_pIC50_direct_inhibition",
    "CYP2C9_pIC50_direct_inhibition",
    "CYP3A4_pIC50_direct_inhibition",
]
for endpoint in endpoints:
    short = endpoint.split("_")[0].lower()
    config = root / f"configs/regression/{short}-molformer-target-transfer-v1.json"
    payload = {
        "track": "regression",
        "endpoint": endpoint,
        "split_run": "20260911T104703655337Z-a3923255",
        "feature_mode": "molformer_pooling",
        "encoder_run": "20260914T140517644424Z-de551457",
        "feature_views": {
            "scores_molformer_cls_mean": "scores47 + CLS768 + mean768 + nested log2fc4"
        },
        "methods": ["random_forest", "extra_trees_l1", "catboost_mae"],
        "blend_weights": [0.0, 0.1, 0.25, 0.5, 1.0],
        "tree_estimators": 240,
        "threads": 4,
        "seed": 20260914,
        "evaluation": "D6-selected view and fixed tree recipes transferred without retuning.",
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
        f"Transfer the fixed MoLFormer tree representation to {endpoint}",
        "--config",
        str(config),
    ]
    for path in inputs:
        command += ["--input", str(path)]
    subprocess.run(
        command + ["--", "python", "scripts/regression_d6_generic_encoder_screen.py"],
        cwd=root,
        check=True,
    )
