"""Record the seven accepted classification states and their LB changes."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
config = root / "configs/analysis/classification-historical-lb-inventory-v1.json"
submission_root = root / "artifacts/submissions/classification"
api_root = root / "artifacts/api/classification"
submission_ids = {
    "0910": "20260910T133122668124Z-68969a86",
    "0912": "20260912T122715261670Z-363bf7b8",
    "0913": "20260913T113807374506Z-77f4825a",
    "0914": "20260914T020759952207Z-37dd7bd8",
    "0915": "20260915T130004513006Z-75b25070",
    "0918": "20260918T094451131332Z-f23dd61b",
    "0923": "20260923T180233070503Z-552d9af4",
}
leaderboard_ids = {
    "0910": ("20260910T135018498829Z-6408f824", "20260910T135018507465Z-d8cffd06"),
    "0912": ("20260912T124752635013Z-591d549f", "20260912T124752634908Z-3ddcd27d"),
    "0913": ("20260913T115952558598Z-27dbd92d", "20260913T115956532758Z-d04c1f66"),
    "0914": ("20260914T021850504300Z-f8bc9ee4", "20260914T021854541658Z-77f4b6cd"),
    "0915": ("20260915T130816907483Z-fa8e6053", "20260915T130816924118Z-19f008d0"),
    "0918": ("20260918T101216302522Z-4a8b944c", "20260918T101216302525Z-527127ee"),
    "0923": ("20260923T183846404373Z-b2478be3", "20260923T183904488705Z-28728897"),
}
inputs = {}
for date, submission_id in submission_ids.items():
    inputs[f"submission-{date}.parquet"] = submission_root / submission_id / "predictions.parquet"
    for endpoint, leaderboard_id in zip(("d6", "a4"), leaderboard_ids[date], strict=True):
        matches = list((api_root / leaderboard_id / "downloads").rglob("*.csv"))
        if len(matches) != 1:
            raise ValueError(f"Expected one leaderboard CSV for {leaderboard_id}")
        inputs[f"{endpoint}-{date}.csv"] = matches[0]

declared = json.loads(config.read_text())
required = {
    state[key]
    for state in declared["states"]
    for key in ("submission", "d6_lb", "a4_lb")
}
if set(inputs) != required:
    raise ValueError("Input aliases differ from config")
for path in inputs.values():
    if not path.is_file():
        raise FileNotFoundError(path)
command = [
    "pixi", "run", "cyp", "run", "--track", "classification",
    "--purpose", "Freeze exact accepted TDI vectors and historical endpoint LB movements",
    "--config", str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_historical_lb_inventory.py"],
    cwd=root,
    check=True,
)
