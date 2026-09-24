"""Launch the no-manual-weight TabPFN 3.5 D6 bank extension."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
inputs = {
    "legacy-bank.parquet": runs
    / "20260915T121009912691Z-73198c15/outputs/prediction-bank.parquet",
    "deployed-oof.parquet": runs
    / "20260918T063800192024Z-c1058115/outputs/oof.parquet",
    "tabpfn35-fast-oof.parquet": runs
    / "20260918T063048677015Z-e0583fe4/outputs/oof.parquet",
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
    "Add TabPFN 3.5-fast to the CYP2D6 prediction bank without manual blend weights",
    "--config",
    str(root / "configs/classification/d6-tabpfn35-bank-extension-v1.json"),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_tabpfn35_bank_extension.py"],
    cwd=root,
    check=True,
)
