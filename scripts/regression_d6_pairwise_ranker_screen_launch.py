"""Register and launch the CYP2D6 pairwise-ranking model screen."""

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

config = root / "configs/regression/d6-pairwise-ranker-screen-v1.json"
payload = {
    "track": "regression",
    "endpoint": "CYP2D6_pIC50_direct_inhibition",
    "split_run": "20260911T104703655337Z-a3923255",
    "views": ["scores51", "medium604", "hybrid2699"],
    "methods": [
        "pairwise_logistic",
        "rank_svm",
        "ranknet",
        "xgb_pairwise",
        "cat_pairlogit",
    ],
    "compact_methods": ["pairwise_logistic", "rank_svm", "ranknet"],
    "maximum_pairs": 30000,
    "pca_dimensions": 64,
    "blend_weights": [0.0, 0.1, 0.25, 0.5, 1.0],
    "threads": 4,
    "ranknet": {
        "epochs": 60,
        "batch_size": 512,
        "learning_rate": 0.001,
        "weight_decay": 0.001,
    },
    "evaluation": (
        "Fixed parent CV; interval-confident pair construction inside outer training only; "
        "positive training-only affine calibration; fixed blends with latest development OOF."
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
        "Compare interval-confident linear, SVM, RankNet, XGBoost, and CatBoost pairwise "
        "CYP2D6 rankers under retained parent CV"
    ),
    "--config",
    str(config),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/regression_d6_pairwise_ranker_screen.py"],
    cwd=root,
    check=True,
)
