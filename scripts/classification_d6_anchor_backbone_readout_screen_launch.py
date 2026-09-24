"""Launch heterogeneous readouts for the D6 anchor/CheMeleon backbone."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
views = runs / "regression/20260917T111518065286Z-801ab060/outputs"
config = root / "configs/classification/d6-anchor-backbone-readout-screen-v1.json"
inputs = {
    "measurements.parquet": runs
    / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    "test-mapping.parquet": runs
    / "shared/20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
    "assignments.parquet": runs
    / "shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    "jeremy-raw.parquet": runs
    / "shared/20260912T055819206227Z-c060d197/outputs/external-embeddings.parquet",
    "jeremy-scores.parquet": runs
    / "shared/20260912T072105522912Z-53eb675c/outputs/external-scores.parquet",
    "chemeleon-raw.parquet": runs
    / "shared/20260912T061231617816Z-299439e8/outputs/external-embeddings.parquet",
    "anchor-raw-scores.parquet": views / "jeremy_init_task_anchor-raw-scores.parquet",
    "pub030-raw.parquet": views / "jeremy_pubchem_anchor030-raw.parquet",
}
declared = json.loads(config.read_text())
needed = {spec["file"] for spec in declared["blocks"].values()}
needed.update({declared["labels"], declared["test_mapping"], declared["assignments"]})
if needed != set(inputs):
    raise ValueError(f"Input mismatch: missing={needed-set(inputs)}, extra={set(inputs)-needed}")
for path in inputs.values():
    if not path.is_file():
        raise FileNotFoundError(path)

command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    (
        "Build heterogeneous classifiers around the new anchored Jeremy and "
        "CheMeleon CYP2D6 backbone"
    ),
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_anchor_backbone_readout_screen.py"],
    cwd=root,
    check=True,
)
