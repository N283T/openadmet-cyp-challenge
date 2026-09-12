"""Package the target-wise regression candidate with top-2 CYP3A4 Caruana."""

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
    "measurements": add(
        "artifacts/runs/regression/20260912T074740645032Z-45d23980/inputs/000-measurements.parquet"
    ),
    "aliases": add(
        "artifacts/runs/regression/20260912T074740645032Z-45d23980/inputs/001-aliases.parquet"
    ),
    "reference": add(
        "artifacts/runs/regression/20260912T074740645032Z-45d23980/inputs/002-cyp-challenge-TEST-BLINDED.csv"
    ),
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
        "predictions": add(f"artifacts/runs/regression/{run}/outputs/test-predictions.parquet"),
        "membership": add(f"artifacts/runs/regression/{run}/outputs/membership.parquet"),
    }
    for name, run in component_runs
]
cfg.update(
    endpoint_components={
        "CYP1A2_pIC50_direct_inhibition": [
            "caruana_3a4_top2",
            "jeremy_frozen_tabpfn",
            "jeremy_scores_log2fc_tabpfn",
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
        "CYP3A4_pIC50_direct_inhibition": [
            "caruana_3a4_top2",
            "jeremy_frozen_tabpfn",
            "jeremy_scores_log2fc_tabpfn",
            "jeremy_finetuned_log2fc_lgbm",
            "chemeleon_log2fc_tabicl",
        ],
    },
    selection_run="20260912T150626214420Z-d602ffed",
    selection_rule=(
        "Use the prior target-wise top-level component selection; within the CYP3A4 "
        "Caruana component use the two highest aggregate-weight full-data models, "
        "renormalized to sum to one."
    ),
    internal_model_count_before=77,
    internal_model_count_after=59,
    caruana_pruning_evidence_run="20260912T152832952478Z-3e4de7c8",
)

config = root / "configs/regression/component-pruned-caruana-3a4-top2-package-v1.json"
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
    (
        "Package the target-wise regression candidate with CYP3A4 Caruana "
        "compressed from 20 models to two"
    ),
    "--config",
    str(config),
]
for path in inputs:
    command.extend(["--input", str(path)])
subprocess.run(
    command + ["--", "python", "scripts/package_pruned_regression.py"],
    check=True,
    cwd=root,
)
