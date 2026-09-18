"""Add each target-leading feature-combination LGBM to its current flat bank."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/regression"
atomic = runs / "20260916T065043779157Z-99517b89/outputs"
broad = runs / "20260917T111705603506Z-e23f10b5/outputs"
diverse = runs / "20260917T112238487125Z-7178add0/outputs"
fast = runs / "20260917T121017505228Z-433d9acb/outputs"
full = runs / "20260917T122531642138Z-54f25a12/outputs"
tdi = runs / "20260917T130425202627Z-93d39ed7/outputs"
feature_screen = runs / "20260918T012518053809Z-5c46e72d/outputs"

inputs = {
    "jeremy-frozen-oof.parquet": runs
    / "20260912T060059705105Z-2f0fbfbd/outputs/external-oof.parquet",
    "jeremy-scores-oof.parquet": runs
    / "20260912T072303592078Z-07f2749a/outputs/external-oof.parquet",
    "jeremy-finetuned-oof.parquet": runs
    / "20260912T072056445264Z-e6b6cbcb/outputs/external-oof.parquet",
    "chemeleon-pca-oof.parquet": runs
    / "20260912T071619253798Z-73a6fb7e/outputs/external-oof.parquet",
    "chemeleon-raw-oof.parquet": runs
    / "20260916T060715739725Z-b0ddf3c8/outputs/external-oof.parquet",
    "owned-2c9-oof.parquet": runs
    / "20260916T095256726141Z-5ba314c0/outputs/external-oof.parquet",
    "3a4-class-quantile-oof.parquet": atomic / "atomic-class_quantile_0.65.parquet",
    "3a4-morgan-quantile-oof.parquet": atomic / "atomic-morgan_quantile_0.65.parquet",
    "owned-raw-scores-lgbm.parquet": broad / "owned_init_hybrid_raw_scores__lgbm-oof.parquet",
    "anchor59-raw-scores-lgbm.parquet": broad
    / "jeremy_init_task_anchor_raw_scores__lgbm-oof.parquet",
    "anchor59-scores-lgbm.parquet": broad
    / "jeremy_init_task_anchor_scores__lgbm-oof.parquet",
    "pubonly-raw-lgbm.parquet": broad / "jeremy_pubchem_only_raw__lgbm-oof.parquet",
    "scratch-raw-catboost.parquet": diverse
    / "scratch_teacher_only_raw__catboost-oof.parquet",
    "anchor59-raw-scores-catboost.parquet": diverse
    / "jeremy_init_task_anchor_raw_scores__catboost-oof.parquet",
    "t35-fast-pub030-raw.parquet": fast / "pub030_raw__tabpfn-oof.parquet",
    "t35-fast-jeremy-raw.parquet": fast / "jeremy_raw__tabpfn-oof.parquet",
    "t35-fast-jeremy-scores-logfc.parquet": fast
    / "jeremy_scores_logfc__tabpfn-oof.parquet",
    "t35-fast-anchor59-raw-scores-logfc.parquet": fast
    / "anchor59_raw_scores_logfc__tabpfn-oof.parquet",
    "t35-fast-chemeleon-raw-logfc.parquet": fast
    / "chemeleon_raw_logfc__tabpfn-oof.parquet",
    "t35-full-pub030-raw.parquet": full / "pub030_raw__tabpfn-oof.parquet",
    "t35-full-jeremy-raw.parquet": full / "jeremy_raw__tabpfn-oof.parquet",
    "t35-full-jeremy-scores-logfc.parquet": full
    / "jeremy_scores_logfc__tabpfn-oof.parquet",
    "t35-full-anchor59-raw-scores-logfc.parquet": full
    / "anchor59_raw_scores_logfc__tabpfn-oof.parquet",
    "t35-full-chemeleon-raw-logfc.parquet": full
    / "chemeleon_raw_logfc__tabpfn-oof.parquet",
    "tdi-extratrees-oof.parquet": tdi / "member_extratrees-oof.parquet",
    "tdi-lgbm-oof.parquet": tdi / "member_lgbm-oof.parquet",
    "tdi-catboost-oof.parquet": tdi / "member_catboost-oof.parquet",
    "tdi-continuous-arm-oof.parquet": runs
    / "20260912T093415808252Z-708fd656/outputs/external-oof.parquet",
    "tdi-continuous-delta-oof.parquet": runs
    / "20260912T093615055236Z-2a25a497/outputs/external-oof.parquet",
    "lgbm-all-pool-top512.parquet": feature_screen / "all_pool_top512-oof.parquet",
    "lgbm-owned-scratch-anchor.parquet": feature_screen
    / "owned_scratch_anchor-oof.parquet",
    "lgbm-jeremy-anchor-all-log2fc.parquet": feature_screen
    / "jeremy_anchor_all_log2fc-oof.parquet",
    "lgbm-dense-pool-top256.parquet": feature_screen / "dense_pool_top256-oof.parquet",
    "current-control.parquet": runs
    / "20260916T062147025088Z-1d1208a1/outputs/aligned-oof.parquet",
    "custom_scoring_functions.py": root / "artifacts/metric-reference/custom_scoring_functions.py",
}

for target in ("1a2", "2c9", "2d6", "3a4"):
    config = root / f"configs/regression/flat-bank-{target}-lgbm-features-v6.json"
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        (
            "Measure the incremental value of the target-leading LGBM feature "
            f"combination in the {target} flat bank"
        ),
        "--config",
        str(config),
    ]
    for alias, path in inputs.items():
        if path.is_file():
            command.extend(["--input", f"{alias}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/regression_flat_ensemble_bank.py"],
        cwd=root,
        check=True,
    )
