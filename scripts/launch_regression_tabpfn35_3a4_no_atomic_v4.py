"""Audit the CYP3A4 flat bank without legacy Caruana-derived atomic models."""

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "artifacts/runs/regression"
FAST = RUNS / "20260917T121017505228Z-433d9acb/outputs"
FULL = RUNS / "20260917T122531642138Z-54f25a12/outputs"
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
    "t35-fast-anchor59-raw-scores-logfc.parquet": FAST
    / "anchor59_raw_scores_logfc__tabpfn-oof.parquet",
    "t35-full-anchor59-raw-scores-logfc.parquet": FULL
    / "anchor59_raw_scores_logfc__tabpfn-oof.parquet",
    "t35-fast-chemeleon-raw-logfc.parquet": FAST / "chemeleon_raw_logfc__tabpfn-oof.parquet",
    "t35-full-chemeleon-raw-logfc.parquet": FULL / "chemeleon_raw_logfc__tabpfn-oof.parquet",
    "t35-fast-jeremy-scores-logfc.parquet": FAST / "jeremy_scores_logfc__tabpfn-oof.parquet",
    "current-control.parquet": RUNS / "20260916T062147025088Z-1d1208a1/outputs/aligned-oof.parquet",
    "custom_scoring_functions.py": ROOT / "artifacts/metric-reference/custom_scoring_functions.py",
}
config = ROOT / "configs/regression/flat-bank-3a4-tabpfn35-no-atomic-v4.json"
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Audit CYP3A4 TabPFN 3.5 flat bank without legacy Caruana atomics",
    "--config",
    str(config),
]
for alias in sorted(inputs):
    command.extend(["--input", f"{alias}={inputs[alias]}"])
subprocess.run(
    command + ["--", "python", "scripts/regression_flat_ensemble_bank.py"], cwd=ROOT, check=True
)
