"""Append cross-fitted TDI readouts to the CYP3A4 flat regression bank."""

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "artifacts/runs/regression"
ATOMIC = RUNS / "20260916T065043779157Z-99517b89/outputs"
FAST = RUNS / "20260917T121017505228Z-433d9acb/outputs"
FULL = RUNS / "20260917T122531642138Z-54f25a12/outputs"
TDI = RUNS / "20260917T130425202627Z-93d39ed7/outputs"
inputs = {
    "jeremy-frozen-oof.parquet": RUNS
    / "20260912T060059705105Z-2f0fbfbd/outputs/external-oof.parquet",
    "jeremy-scores-oof.parquet": RUNS
    / "20260912T072303592078Z-07f2749a/outputs/external-oof.parquet",
    "jeremy-finetuned-oof.parquet": RUNS
    / "20260912T072056445264Z-e6b6cbcb/outputs/external-oof.parquet",
    "chemeleon-pca-oof.parquet": RUNS
    / "20260912T071619253798Z-73a6fb7e/outputs/external-oof.parquet",
    "chemeleon-raw-oof.parquet": RUNS
    / "20260916T060715739725Z-b0ddf3c8/outputs/external-oof.parquet",
    "3a4-class-quantile-oof.parquet": ATOMIC / "atomic-class_quantile_0.65.parquet",
    "3a4-morgan-quantile-oof.parquet": ATOMIC / "atomic-morgan_quantile_0.65.parquet",
    "t35-fast-anchor59-raw-scores-logfc.parquet": FAST
    / "anchor59_raw_scores_logfc__tabpfn-oof.parquet",
    "t35-full-anchor59-raw-scores-logfc.parquet": FULL
    / "anchor59_raw_scores_logfc__tabpfn-oof.parquet",
    "t35-fast-chemeleon-raw-logfc.parquet": FAST / "chemeleon_raw_logfc__tabpfn-oof.parquet",
    "t35-full-chemeleon-raw-logfc.parquet": FULL / "chemeleon_raw_logfc__tabpfn-oof.parquet",
    "t35-fast-jeremy-scores-logfc.parquet": FAST / "jeremy_scores_logfc__tabpfn-oof.parquet",
    "tdi-extratrees-oof.parquet": TDI / "member_extratrees-oof.parquet",
    "tdi-lgbm-oof.parquet": TDI / "member_lgbm-oof.parquet",
    "tdi-catboost-oof.parquet": TDI / "member_catboost-oof.parquet",
    "tdi-continuous-arm-oof.parquet": RUNS
    / "20260912T093415808252Z-708fd656/outputs/external-oof.parquet",
    "tdi-continuous-delta-oof.parquet": RUNS
    / "20260912T093615055236Z-2a25a497/outputs/external-oof.parquet",
    "current-control.parquet": RUNS / "20260916T062147025088Z-1d1208a1/outputs/aligned-oof.parquet",
    "custom_scoring_functions.py": ROOT / "artifacts/metric-reference/custom_scoring_functions.py",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Append TDI probability and continuous-task readouts to the CYP3A4 flat bank",
    "--config",
    str(ROOT / "configs/regression/flat-bank-3a4-tabpfn35-tdi-v5.json"),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_flat_ensemble_bank.py"],
    cwd=ROOT,
    check=True,
)
