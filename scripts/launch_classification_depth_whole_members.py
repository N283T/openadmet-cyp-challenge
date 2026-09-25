"""Record exact submitted TDI OOF reconstruction and whole-depth screen."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
d6 = runs / "20260925T130752164378Z-d40a0d9c"
a4 = runs / "20260915T125001009276Z-d5adff22"
sources = {
    "depth-oof.parquet": root
    / "artifacts/runs/shared/20260925T065626160287Z-8017d516/outputs/oof.parquet",
    "d6-current.parquet": d6 / "outputs/crossfit-mean-block-oof.parquet",
    "d6-lgbm.parquet": d6 / "inputs/lgbm-oof.parquet",
    "d6-pfn.parquet": d6 / "inputs/pfn-oof.parquet",
    "d6-graph.parquet": d6 / "inputs/graph-oof.parquet",
    "d6-fixed.parquet": d6 / "inputs/fixed-oof.parquet",
    "d6-raw600.parquet": d6 / "inputs/raw600-oof.parquet",
    "d6-selections.parquet": d6 / "inputs/selections.parquet",
    "a4-oof.parquet": a4 / "outputs/oof.parquet",
    "a4-member.parquet": a4 / "inputs/000-member-oof.parquet",
    "a4-raw.parquet": a4 / "inputs/004-outer-oof.parquet",
    "a4-aux.parquet": a4 / "inputs/002-oof.parquet",
    "a4-selections.parquet": a4 / "outputs/selections.parquet",
    "tabpfn-oof.parquet": (
        runs / "20260925T140413952333Z-933f232c/outputs/oof-probabilities.parquet"
    ),
    "tabicl-oof.parquet": (
        runs / "20260925T140649186742Z-c7a759af/outputs/outer-oof.parquet"
    ),
}
command = [
    "pixi", "run", "cyp", "run", "--track", "classification", "--purpose",
    "Reconstruct submitted D6 and A4 TDI OOF and screen whole CheMeleon depth members",
    "--config", str(root / "configs/screen/classification-depth-whole-members-v2.json"),
]
for alias, path in sources.items():
    if not path.is_file():
        raise FileNotFoundError(path)
    command.extend(["--input", f"{alias}={path}"])
command.extend(["--", "python", "scripts/evaluate_classification_depth_whole_members.py"])
subprocess.run(command, cwd=root, check=True)
