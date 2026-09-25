"""Record member-level depth swaps against submitted non-D6 regression banks."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/regression"
banks = {
    "1a2": (
        runs / "20260918T034159381060Z-46314b50",
        "outputs/1A2/add_dense_and_fast_pool",
    ),
    "2c9": (runs / "20260917T123601175420Z-e1560701", "outputs"),
    "3a4": (runs / "20260917T130718808881Z-dfa02730", "outputs"),
}
jeremy_runs = {
    "1a2": "20260925T134423571457Z-22c133cf",
    "2c9": "20260925T134513708561Z-b72b5dc0",
    "3a4": "20260925T134540256996Z-d9bc0461",
}
sources = {
    "custom_scoring_functions.py": root / "artifacts/metric-reference/custom_scoring_functions.py",
    "chemeleon-depth-oof.parquet": runs / "20260925T130401763912Z-dcf9b1b8/outputs/oof.parquet",
    "chemeleon-depth1-tabpfn.parquet": (
        runs / "20260925T134457782507Z-c352d23e/outputs/external-oof.parquet"
    ),
    "chemeleon-depth1-tabicl.parquet": (
        runs / "20260925T134743468427Z-a53c9eef/outputs/external-oof.parquet"
    ),
}
for target, (run, output_directory) in banks.items():
    output = run / output_directory
    summary = output / "summary.json"
    sources[f"{target}-summary.json"] = summary
    sources[f"{target}-bank-oof.parquet"] = output / "oof.parquet"
    sources[f"{target}-bank-fits.csv"] = output / "fits.csv"
    sources[f"{target}-jeremy-depth-oof.parquet"] = (
        runs / jeremy_runs[target] / "outputs/oof.parquet"
    )
    for component in json.loads(summary.read_text())["config"]["components"]:
        sources[f"{target}-{component['name']}.parquet"] = run / "inputs" / component["file"]
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Reconstruct submitted 1A2, 2C9 and 3A4 banks and screen depth member swaps",
    "--config",
    str(root / "configs/screen/other-regression-depth-members-v1.json"),
]
for alias, path in sources.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
command.extend(["--", "python", "scripts/evaluate_other_regression_depth_members.py"])
subprocess.run(command, cwd=root, check=True)
