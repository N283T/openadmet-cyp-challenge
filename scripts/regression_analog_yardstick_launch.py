"""Launch the analog-proxy yardstick analysis of the sent regression blends."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
screen = runs / "regression/20261029T145728905917Z-9d7ea648"
inputs = {
    f"control-{name}.parquet": screen / f"inputs/control-{name}.parquet"
    for name in ("1a2", "2c9", "2d6", "3a4")
}
inputs["cypmol-oof.parquet"] = screen / "outputs/predictions/oof.parquet"
inputs["test-mapping.parquet"] = (
    runs / "shared/20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet"
)
inputs["analog-pairs.parquet"] = (
    root / "artifacts/tracked/d6-analog-expansion-validation-20260924/anchor-analog-pairs.parquet"
)
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    (
        "Score the sent regression blends on analog-proxy OOF rows against their Live deltas, "
        "repeat nested selection on the analog objective and scan global shifts"
    ),
    "--config",
    str(root / "configs/regression/analog-yardstick-v1.json"),
]
for alias, path in inputs.items():
    command += ["--input", f"{alias}={path}"]
subprocess.run(
    command + ["--", "python", "scripts/regression_analog_yardstick.py"], cwd=root, check=True
)
