"""Launch fold-safe Jeremy-anchored D6 fine-tuning."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
distillation = runs / "regression/20260917T051118981633Z-a9758b6c/inputs"
inputs = {
    "jeremy-checkpoint.pt": distillation / "jeremy-checkpoint.pt",
    "compounds.parquet": distillation / "compounds.parquet",
    "measurements.parquet": runs
    / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    "assignments.parquet": runs
    / "shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    "test-mapping.parquet": runs
    / "shared/20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
    "chemeleon-raw.parquet": runs
    / "shared/20260912T061231617816Z-299439e8/outputs/external-embeddings.parquet",
}
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
        "Fine-tune Jeremy on fold-safe D6 labels with teacher anchoring and evaluate "
        "adapted raw600 plus frozen CheMeleon readouts"
    ),
    "--config",
    str(root / "configs/classification/d6-jeremy-anchored-finetune-v1.json"),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_jeremy_anchored_finetune.py"],
    cwd=root,
    check=True,
)
