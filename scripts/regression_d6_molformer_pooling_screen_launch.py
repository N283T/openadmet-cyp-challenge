"""Register and launch the MoLFormer pooling/readout CYP2D6 screen."""

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

config = root / "configs/regression/d6-molformer-pooling-screen-v1.json"
payload = {
    "track": "regression",
    "endpoint": "CYP2D6_pIC50_direct_inhibition",
    "split_run": "20260911T104703655337Z-a3923255",
    "feature_mode": "molformer_pooling",
    "encoder_run": "20260914T140517644424Z-de551457",
    "feature_views": {
        "molformer_mean": "masked mean768 + nested log2fc4",
        "molformer_cls_mean": "CLS768 + masked mean768 + nested log2fc4",
        "scores_molformer_mean": "scores47 + masked mean768 + nested log2fc4",
        "scores_molformer_cls_mean": "scores47 + CLS768 + masked mean768 + nested log2fc4",
    },
    "methods": [
        "extra_trees_l1",
        "extra_trees_l3",
        "extra_trees_l8",
        "extra_trees_f025",
        "random_forest",
        "catboost_mae",
        "xgboost",
    ],
    "blend_weights": [0.0, 0.1, 0.25, 0.5, 1.0],
    "tree_estimators": 240,
    "threads": 4,
    "seed": 20260914,
    "evaluation": "Fixed parent CV; fixed pooling and small tree recipes; no tuning on OOF.",
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
    "Compare MoLFormer pooling and regularized tree readouts for CYP2D6 ordering",
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
