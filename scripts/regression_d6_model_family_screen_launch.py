"""Register and launch the broad CYP2D6 fixed model-family screen."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
inputs = [
    runs / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
    runs / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    runs / "shared/20260912T055819206227Z-c060d197/outputs/external-embeddings.parquet",
    runs / "shared/20260912T072105522912Z-53eb675c/outputs/external-scores.parquet",
    runs / "shared/20260912T061043373484Z-2e344d7a/outputs/fp-mordred-features.parquet",
    runs / "shared/20260912T082554377951Z-d272774c/outputs/nested-log2fc.parquet",
    runs / "regression/20260912T083627802028Z-ae562a59/outputs/focus_2d6_combined_50pct.parquet",
    runs / "regression/20260912T151747998390Z-f24622ac/outputs/crossfit-predictions.parquet",
]
for path in inputs:
    if not path.is_file():
        raise FileNotFoundError(path)

config = root / "configs/regression/d6-model-family-screen-v1.json"
payload = {
    "track": "regression",
    "endpoint": "CYP2D6_pIC50_direct_inhibition",
    "split_run": "20260911T104703655337Z-a3923255",
    "feature_views": {
        "scores51": "upstream predictor scores47 + new nested log2fc4",
        "medium604": "frozen upstream medium embedding600 + new nested log2fc4",
        "hybrid2699": "medium600 + scores47 + Morgan2048 + new nested log2fc4",
    },
    "methods": [
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
    ],
    "blend_weights": [0.0, 0.1, 0.25, 0.5, 1.0],
    "tree_estimators": 240,
    "pca_dimensions": 64,
    "threads": 4,
    "kan": {
        "pca_dimensions": 64,
        "hidden": 16,
        "grid_size": 8,
        "bound": 3.0,
        "epochs": 180,
        "batch_size": 128,
        "learning_rate": 0.002,
        "weight_decay": 0.0001,
    },
    "evaluation": (
        "Fixed parent CV; inspect standalone ordering, residual diversity, and prespecified "
        "raw/moment-matched blends with the latest development OOF."
    ),
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
    (
        "Explore classical, modern tabular, kernel, neural, and compact KAN CYP2D6 "
        "readouts for complementary molecular ordering"
    ),
    "--config",
    str(config),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/regression_d6_model_family_screen.py"],
    cwd=root,
    check=True,
)
