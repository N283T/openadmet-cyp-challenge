"""Launch target-specific flat banks for atomic distilled-checkpoint readouts."""

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "artifacts/runs/regression"
BROAD = RUNS / "20260917T111705603506Z-e23f10b5/outputs"
DIVERSE = RUNS / "20260917T112238487125Z-7178add0/outputs"
inputs = {
    "jeremy-medium-oof.parquet": RUNS
    / "20260912T055849013223Z-f16a5dfd/outputs/external-oof.parquet",
    "current-control.parquet": RUNS
    / "20260916T062147025088Z-1d1208a1/outputs/aligned-oof.parquet",
    "custom_scoring_functions.py": ROOT
    / "artifacts/metric-reference/custom_scoring_functions.py",
}
broad = {
    "scratch-raw-lgbm": "scratch_teacher_only_raw__lgbm-oof.parquet",
    "owned-raw-lgbm": "owned_init_hybrid_raw__lgbm-oof.parquet",
    "owned-hidden-lgbm": "owned_init_hybrid_hidden__lgbm-oof.parquet",
    "owned-raw-scores-lgbm": "owned_init_hybrid_raw_scores__lgbm-oof.parquet",
    "anchor59-raw-lgbm": "jeremy_init_task_anchor_raw__lgbm-oof.parquet",
    "anchor59-hidden-lgbm": "jeremy_init_task_anchor_hidden__lgbm-oof.parquet",
    "anchor59-scores-lgbm": "jeremy_init_task_anchor_scores__lgbm-oof.parquet",
    "anchor59-raw-scores-lgbm": "jeremy_init_task_anchor_raw_scores__lgbm-oof.parquet",
    "pubonly-raw-lgbm": "jeremy_pubchem_only_raw__lgbm-oof.parquet",
    "pub010-raw-lgbm": "jeremy_pubchem_anchor010_raw__lgbm-oof.parquet",
    "pub010-raw-scores-lgbm": "jeremy_pubchem_anchor010_raw_scores__lgbm-oof.parquet",
    "pub030-raw-lgbm": "jeremy_pubchem_anchor030_raw__lgbm-oof.parquet",
    "pub030-raw-scores-lgbm": "jeremy_pubchem_anchor030_raw_scores__lgbm-oof.parquet",
}
diverse = {
    "scratch-raw-catboost": "scratch_teacher_only_raw__catboost-oof.parquet",
    "scratch-raw-tabpfn": "scratch_teacher_only_raw__tabpfn-oof.parquet",
    "owned-raw-scores-catboost": "owned_init_hybrid_raw_scores__catboost-oof.parquet",
    "owned-raw-scores-tabpfn": "owned_init_hybrid_raw_scores__tabpfn-oof.parquet",
    "anchor59-raw-scores-catboost": "jeremy_init_task_anchor_raw_scores__catboost-oof.parquet",
    "anchor59-raw-scores-tabpfn": "jeremy_init_task_anchor_raw_scores__tabpfn-oof.parquet",
    "anchor59-scores-tabpfn": "jeremy_init_task_anchor_scores__tabpfn-oof.parquet",
    "pub030-raw-tabpfn": "jeremy_pubchem_anchor030_raw__tabpfn-oof.parquet",
    "pub030-raw-scores-catboost": "jeremy_pubchem_anchor030_raw_scores__catboost-oof.parquet",
}
for alias, filename in broad.items():
    inputs[f"broad-{alias}.parquet"] = BROAD / filename
for alias, filename in diverse.items():
    inputs[f"diverse-{alias}.parquet"] = DIVERSE / filename

for target in ["1a2", "2c9", "2d6", "3a4"]:
    config = ROOT / f"configs/regression/distilled-flat-bank-{target}-v1.json"
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        f"Evaluate target-specific flat ensemble methods for distilled checkpoint bank {target}",
        "--config",
        str(config),
    ]
    for alias, path in inputs.items():
        command.extend(["--input", f"{alias}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/regression_flat_ensemble_bank.py"],
        cwd=ROOT,
        check=True,
    )
