"""Record D6 member addition and replacement comparisons."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/regression"
sources = {
    "current-tabpfn.parquet": (
        runs / "20260912T072303592078Z-07f2749a/outputs/external-oof.parquet"
    ),
    "current-ft-lgbm.parquet": (
        runs / "20260912T072056445264Z-e6b6cbcb/outputs/external-oof.parquet"
    ),
    "chemeleon-depth-oof.parquet": runs / "20260925T130401763912Z-dcf9b1b8/outputs/oof.parquet",
    "chemeleon-d1-tabpfn.parquet": (
        runs / "20260925T132723475153Z-14683186/outputs/external-oof.parquet"
    ),
    "chemeleon-d1-tabicl.parquet": (
        runs / "20260925T132758390039Z-67f7eeec/outputs/external-oof.parquet"
    ),
    "chemeleon-d1-tabpfn-log2fc.parquet": (
        runs / "20260925T133017168391Z-0fdef233/outputs/external-oof.parquet"
    ),
    "chemeleon-d1-tabicl-log2fc.parquet": (
        runs / "20260925T133115393813Z-d56c7302/outputs/external-oof.parquet"
    ),
    "jeremy-depth-oof.parquet": runs / "20260925T133311172709Z-39e028eb/outputs/oof.parquet",
    "custom_scoring_functions.py": root / "artifacts/metric-reference/custom_scoring_functions.py",
}
command = [
    "pixi", "run", "cyp", "run", "--track", "regression",
    "--purpose", "Compare whole-model CYP2D6 depth additions and FT replacement",
    "--config", str(root / "configs/screen/d6-depth-atomic-ensemble-v1.json"),
]
for alias, path in sources.items():
    command.extend(["--input", f"{alias}={path}"])
command.extend(["--", "python", "scripts/evaluate_d6_depth_atomic_ensemble.py"])
subprocess.run(command, cwd=root, check=True)
