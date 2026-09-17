"""Append checkpoint-derived atomic predictions to the established CYP2D6 bank."""

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "artifacts/runs/regression"
BROAD = RUNS / "20260917T111705603506Z-e23f10b5/outputs"
DIVERSE = RUNS / "20260917T112238487125Z-7178add0/outputs"
inputs = {
    "jeremy-frozen-oof.parquet": RUNS
    / "20260912T060059705105Z-2f0fbfbd/outputs/external-oof.parquet",
    "jeremy-scores-oof.parquet": RUNS
    / "20260912T072303592078Z-07f2749a/outputs/external-oof.parquet",
    "jeremy-finetuned-oof.parquet": RUNS
    / "20260912T072056445264Z-e6b6cbcb/outputs/external-oof.parquet",
    "chemeleon-raw-oof.parquet": RUNS
    / "20260916T060715739725Z-b0ddf3c8/outputs/external-oof.parquet",
    "current-control.parquet": RUNS
    / "20260916T062147025088Z-1d1208a1/outputs/aligned-oof.parquet",
    "custom_scoring_functions.py": ROOT
    / "artifacts/metric-reference/custom_scoring_functions.py",
    "owned-raw-scores-lgbm.parquet": BROAD
    / "owned_init_hybrid_raw_scores__lgbm-oof.parquet",
    "anchor59-raw-scores-lgbm.parquet": BROAD
    / "jeremy_init_task_anchor_raw_scores__lgbm-oof.parquet",
    "anchor59-scores-lgbm.parquet": BROAD
    / "jeremy_init_task_anchor_scores__lgbm-oof.parquet",
    "pubonly-raw-lgbm.parquet": BROAD / "jeremy_pubchem_only_raw__lgbm-oof.parquet",
    "scratch-raw-catboost.parquet": DIVERSE
    / "scratch_teacher_only_raw__catboost-oof.parquet",
    "anchor59-raw-scores-catboost.parquet": DIVERSE
    / "jeremy_init_task_anchor_raw_scores__catboost-oof.parquet",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Append distilled checkpoint atomic predictions to the established CYP2D6 flat bank",
    "--config",
    str(ROOT / "configs/regression/flat-bank-2d6-distilled-v2.json"),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_flat_ensemble_bank.py"],
    cwd=ROOT,
    check=True,
)
