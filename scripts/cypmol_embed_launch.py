"""Launch the frozen CYPMol inhibitor-model embedding export in the cypmol pixi environment."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
ghq = Path.home() / "ghq/github.com"
inputs = {
    "compounds.parquet": root
    / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    "dict.txt": ghq / "deepmodeling/Uni-Mol/unimol/example_data/molecule/dict.txt",
    "inhibitor.zip": ghq / "CjmTH/CYPMol/Inhibitor/inhibitor.zip",
    "cypmol_model.py": ghq / "CjmTH/CYPMol/Inhibitor/main_seq_pocket_predict.py",
}
members = (0, 1, 2, 3, 4, 7, 9, 10, 12, 15)
inputs.update(
    {f"checkpoint-sub{n}.pt": Path.home() / f".cache/cypmol/inhi_seq_pocket_sub{n}.pt" for n in members}
)
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "shared",
    "--purpose",
    (
        "Export frozen CYPMol (CjmTH/CYPMol 0c6657d, inhibitor ensemble of 10) 512-d molecule-by-isoform "
        "vectors and inhibitor probabilities for all 6896 compounds, mean over 10 ensemble members and three conformers"
    ),
    "--config",
    str(root / "configs/shared/cypmol-embed-ensemble-v1.json"),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", str(root / ".pixi/envs/cypmol/bin/python"), "scripts/cypmol_embed.py"],
    cwd=root,
    check=True,
)
