"""Append selected TabPFN 3.5 peers and rerun each target-specific flat bank."""

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "artifacts/runs/regression"
ATOMIC_3A4 = RUNS / "20260916T065043779157Z-99517b89/outputs"
BROAD = RUNS / "20260917T111705603506Z-e23f10b5/outputs"
DIVERSE = RUNS / "20260917T112238487125Z-7178add0/outputs"
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
    "owned-2c9-oof.parquet": RUNS / "20260916T095256726141Z-5ba314c0/outputs/external-oof.parquet",
    "3a4-class-quantile-oof.parquet": ATOMIC_3A4 / "atomic-class_quantile_0.65.parquet",
    "3a4-morgan-quantile-oof.parquet": ATOMIC_3A4 / "atomic-morgan_quantile_0.65.parquet",
    "owned-raw-scores-lgbm.parquet": BROAD / "owned_init_hybrid_raw_scores__lgbm-oof.parquet",
    "anchor59-raw-scores-lgbm.parquet": BROAD
    / "jeremy_init_task_anchor_raw_scores__lgbm-oof.parquet",
    "anchor59-scores-lgbm.parquet": BROAD / "jeremy_init_task_anchor_scores__lgbm-oof.parquet",
    "pubonly-raw-lgbm.parquet": BROAD / "jeremy_pubchem_only_raw__lgbm-oof.parquet",
    "scratch-raw-catboost.parquet": DIVERSE / "scratch_teacher_only_raw__catboost-oof.parquet",
    "anchor59-raw-scores-catboost.parquet": DIVERSE
    / "jeremy_init_task_anchor_raw_scores__catboost-oof.parquet",
    "t35-fast-pub030-raw.parquet": FAST / "pub030_raw__tabpfn-oof.parquet",
    "t35-fast-jeremy-raw.parquet": FAST / "jeremy_raw__tabpfn-oof.parquet",
    "t35-fast-jeremy-scores-logfc.parquet": FAST / "jeremy_scores_logfc__tabpfn-oof.parquet",
    "t35-fast-anchor59-raw-scores-logfc.parquet": FAST
    / "anchor59_raw_scores_logfc__tabpfn-oof.parquet",
    "t35-fast-chemeleon-raw-logfc.parquet": FAST / "chemeleon_raw_logfc__tabpfn-oof.parquet",
    "t35-full-pub030-raw.parquet": FULL / "pub030_raw__tabpfn-oof.parquet",
    "t35-full-jeremy-raw.parquet": FULL / "jeremy_raw__tabpfn-oof.parquet",
    "t35-full-jeremy-scores-logfc.parquet": FULL / "jeremy_scores_logfc__tabpfn-oof.parquet",
    "t35-full-anchor59-raw-scores-logfc.parquet": FULL
    / "anchor59_raw_scores_logfc__tabpfn-oof.parquet",
    "t35-full-chemeleon-raw-logfc.parquet": FULL / "chemeleon_raw_logfc__tabpfn-oof.parquet",
    "current-control.parquet": RUNS / "20260916T062147025088Z-1d1208a1/outputs/aligned-oof.parquet",
    "custom_scoring_functions.py": ROOT / "artifacts/metric-reference/custom_scoring_functions.py",
}

for target in ("1a2", "2c9", "2d6", "3a4"):
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        f"Rerun {target} flat atomic bank with TabPFN 3.5 fast/full peers",
        "--config",
        str(ROOT / f"configs/regression/flat-bank-{target}-tabpfn35-v3.json"),
    ]
    for alias, path in inputs.items():
        if path.exists():
            command.extend(["--input", f"{alias}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/regression_flat_ensemble_bank.py"], cwd=ROOT, check=True
    )
