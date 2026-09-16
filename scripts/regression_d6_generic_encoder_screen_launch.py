"""Register and launch the generic-encoder CYP2D6 screen."""

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
    runs / "regression/20260914T135108660328Z-387c0bcc/outputs/features.parquet",
    runs / "regression/20260914T135140332924Z-1d133b02/outputs/features.parquet",
]
for path in inputs:
    if not path.is_file():
        raise FileNotFoundError(path)

config = root / "configs/regression/d6-generic-encoder-screen-v1.json"
payload = {
    "track": "regression",
    "endpoint": "CYP2D6_pIC50_direct_inhibition",
    "split_run": "20260911T104703655337Z-a3923255",
    "encoder_runs": {
        "kermt_v2": "20260914T135108660328Z-387c0bcc",
        "molformer_c3": "20260914T135140332924Z-1d133b02",
    },
    "feature_views": {
        "kermt": "generic KERMT v2 3200 + nested log2fc4",
        "molformer": "generic MoLFormer-c3 768 + nested log2fc4",
        "scores_kermt": "upstream scores47 + generic KERMT v2 3200 + nested log2fc4",
        "scores_molformer": "upstream scores47 + generic MoLFormer-c3 768 + nested log2fc4",
    },
    "methods": ["pca64_ridge", "pca64_rbf_svr", "pls16", "extra_trees"],
    "blend_weights": [0.0, 0.1, 0.25, 0.5, 1.0],
    "tree_estimators": 300,
    "threads": 4,
    "seed": 20260914,
    "evaluation": "Fixed parent CV and fixed readouts; no tuning on this representation screen.",
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
    "Test untouched generic KERMT v2 and MoLFormer-c3 representations for CYP2D6 ordering",
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
