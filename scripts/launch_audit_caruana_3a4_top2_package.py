"""Record the independent audit for the CYP3A4 top-2 package."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
inputs = []


def add(relative: str) -> str:
    path = root / relative
    inputs.append(path)
    return f"{len(inputs) - 1:03d}-{path.name}"


cfg = {
    "packaged_long": add(
        "artifacts/runs/regression/20260912T163036856426Z-bf2e35dc/outputs/test-predictions.parquet"
    ),
    "packaged_wide": add(
        "artifacts/runs/regression/20260912T163036856426Z-bf2e35dc/outputs/regression.parquet"
    ),
    "previous_packaged_long": add(
        "artifacts/runs/regression/20260912T151847989358Z-14af08d1/outputs/test-predictions.parquet"
    ),
    "reference": add(
        "artifacts/runs/regression/20260912T074740645032Z-45d23980/inputs/002-cyp-challenge-TEST-BLINDED.csv"
    ),
    "reconstruction": add(
        "artifacts/runs/regression/20260912T152832952478Z-3e4de7c8/outputs/reconstruction.csv"
    ),
    "pruning_scores": add(
        "artifacts/runs/regression/20260912T152832952478Z-3e4de7c8/outputs/scores.csv"
    ),
    "top_level_scores": add(
        "artifacts/runs/regression/20260912T151747998390Z-f24622ac/outputs/crossfit-scores.csv"
    ),
    "endpoint": "CYP3A4_pIC50_direct_inhibition",
    "unchanged_endpoints": [
        "CYP1A2_pIC50_direct_inhibition",
        "CYP2C9_pIC50_direct_inhibition",
        "CYP2D6_pIC50_direct_inhibition",
    ],
    "top_level_counts": {
        "CYP1A2_pIC50_direct_inhibition": 3,
        "CYP2C9_pIC50_direct_inhibition": 3,
        "CYP2D6_pIC50_direct_inhibition": 2,
        "CYP3A4_pIC50_direct_inhibition": 5,
    },
    "caruana_internal_models_after": 2,
}
component_runs = [
    ("caruana_3a4_top2", "20260912T153211990718Z-e29d7578"),
    ("jeremy_frozen_tabpfn", "20260912T074409483212Z-f3024998"),
    ("jeremy_scores_log2fc_tabpfn", "20260912T074557646349Z-da857069"),
    ("jeremy_finetuned_log2fc_lgbm", "20260912T074624192288Z-bd58b762"),
    ("chemeleon_log2fc_tabicl", "20260912T074558733332Z-98ab8ceb"),
]
cfg["components"] = [
    {
        "name": name,
        "test": add(f"artifacts/runs/regression/{run}/outputs/test-predictions.parquet"),
    }
    for name, run in component_runs
]
cfg["endpoint_components"] = {
    "CYP1A2_pIC50_direct_inhibition": [
        "jeremy_frozen_tabpfn",
        "jeremy_scores_log2fc_tabpfn",
        "jeremy_finetuned_log2fc_lgbm",
    ],
    "CYP2C9_pIC50_direct_inhibition": [
        "jeremy_scores_log2fc_tabpfn",
        "jeremy_finetuned_log2fc_lgbm",
        "chemeleon_log2fc_tabicl",
    ],
    "CYP2D6_pIC50_direct_inhibition": [
        "jeremy_scores_log2fc_tabpfn",
        "jeremy_finetuned_log2fc_lgbm",
    ],
    "CYP3A4_pIC50_direct_inhibition": [name for name, _run in component_runs],
}

config = root / "configs/regression/component-pruned-caruana-3a4-top2-audit-v3.json"
with config.open("x") as stream:
    json.dump(cfg, stream, indent=2)
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Independently audit the target-wise regression package with top-2 CYP3A4 Caruana",
    "--config",
    str(config),
]
for path in inputs:
    command.extend(["--input", str(path)])
subprocess.run(
    command + ["--", "python", "scripts/audit_caruana_3a4_top2_package.py"],
    check=True,
    cwd=root,
)
