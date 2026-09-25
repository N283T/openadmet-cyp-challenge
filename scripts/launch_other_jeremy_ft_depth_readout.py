"""Record FT600 depth refits for one non-D6 regression endpoint."""

import argparse
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("target", choices=["1a2", "2c9", "3a4"])
target = parser.parse_args().target
endpoint = f"CYP{target.upper()}_pIC50_direct_inhibition"
runs = root / "artifacts/runs/regression"
depth = runs / "20260925T132940652082Z-227d97fa/outputs"
original = runs / "20260912T072056445264Z-e6b6cbcb/outputs"
sources = {
    "parent-assignments.parquet": (
        runs / "20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet"
    ),
    "measurements.parquet": (
        root
        / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
        / "measurements.parquet"
    ),
    "nested-log2fc.parquet": (
        root
        / "artifacts/runs/shared/20260912T052549384296Z-bac5fb3a/outputs"
        / "nested-log2fc.parquet"
    ),
    "original-ft-oof.parquet": original / "external-oof.parquet",
    "custom_scoring_functions.py": root / "artifacts/metric-reference/custom_scoring_functions.py",
}
for fold in range(5):
    sources[f"fold{fold}-original-final.parquet"] = (
        runs / "20260912T060538153709Z-8d434eee/outputs" / f"fold{fold}/external-embeddings.parquet"
    )
    sources[f"fold{fold}-selected.json"] = original / f"fold-{fold}/{endpoint}/selected.json"
    for level in (1, 2, 4):
        sources[f"fold{fold}-depth{level}.parquet"] = depth / f"fold{fold}-depth{level}.parquet"
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    (f"CYP{target.upper()} Jeremy FT600 depth refit with unchanged inner-selected LightGBM"),
    "--config",
    str(root / f"configs/screen/jeremy-ft-depth-readout-{target}-v1.json"),
]
for alias, path in sources.items():
    command.extend(["--input", f"{alias}={path}"])
command.extend(["--", "python", "scripts/evaluate_d6_jeremy_ft_depths.py"])
subprocess.run(command, cwd=root, check=True)
