"""Record 2C9/3A4 bank reselection with whole depth members."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/regression"
banks = {
    "2c9": runs / "20260917T123601175420Z-e1560701",
    "3a4": runs / "20260917T130718808881Z-dfa02730",
}
jeremy_runs = {
    "2c9": "20260925T134513708561Z-b72b5dc0",
    "3a4": "20260925T134540256996Z-d9bc0461",
}
sources = {
    "custom_scoring_functions.py": root / "artifacts/metric-reference/custom_scoring_functions.py",
    "chemeleon-depth-oof.parquet": runs / "20260925T130401763912Z-dcf9b1b8/outputs/oof.parquet",
    "chemeleon-depth1-tabicl.parquet": (
        runs / "20260925T134743468427Z-a53c9eef/outputs/external-oof.parquet"
    ),
}
for target, run in banks.items():
    output = run / "outputs"
    summary = output / "summary.json"
    sources[f"{target}-summary.json"] = summary
    sources[f"{target}-bank-oof.parquet"] = output / "oof.parquet"
    sources[f"{target}-jeremy-depth-oof.parquet"] = (
        runs / jeremy_runs[target] / "outputs/oof.parquet"
    )
    for component in json.loads(summary.read_text())["config"]["components"]:
        sources[f"{target}-{component['name']}.parquet"] = run / "inputs" / component["file"]
command = [
    "pixi", "run", "cyp", "run", "--track", "regression", "--purpose",
    "Reselect submitted 2C9 and 3A4 banks with whole intermediate-depth members",
    "--config", str(root / "configs/screen/regression-depth-bank-reselection-v1.json"),
]
for alias, path in sources.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
command.extend(["--", "python", "scripts/evaluate_regression_depth_bank_reselection.py"])
subprocess.run(command, cwd=root, check=True)
