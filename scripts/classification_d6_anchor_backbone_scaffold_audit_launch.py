"""Launch scaffold confirmation of UMAP-prespecified D6 new-family subsets."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
config = root / "configs/classification/d6-anchor-backbone-scaffold-audit-v1.json"
inputs = {
    "lgbm-oof.parquet": runs
    / "20260919T155849130417Z-0394b600/outputs/oof.parquet",
    "readout-oof.parquet": runs
    / "20260923T161449619218Z-a5f9f7ea/outputs/oof.parquet",
    "tabpfn-oof.parquet": runs
    / "20260923T161718040624Z-8a9676e7/outputs/oof.parquet",
}
declared = json.loads(config.read_text())
needed = {declared[key] for key in ["lgbm_oof", "readout_oof", "tabpfn_oof"]}
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
    "Confirm UMAP-prespecified D6 anchor/CheMeleon bank subsets on scaffold folds",
    "--config",
    str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command
    + ["--", "python", "scripts/classification_d6_anchor_backbone_scaffold_audit.py"],
    cwd=root,
    check=True,
)
