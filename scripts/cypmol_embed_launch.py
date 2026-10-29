"""Launch the frozen CYPMol inhibitor-model embedding export in the cypmol pixi environment."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
ghq = Path.home() / "ghq/github.com"
inputs = {
    "compounds.parquet": root
    / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    "checkpoint.pt": Path.home() / ".cache/cypmol/inhi_seq_pocket_sub0.pt",
    "dict.txt": ghq / "deepmodeling/Uni-Mol/unimol/example_data/molecule/dict.txt",
    "inhibitor.zip": ghq / "CjmTH/CYPMol/Inhibitor/inhibitor.zip",
    "cypmol_model.py": ghq / "CjmTH/CYPMol/Inhibitor/main_seq_pocket_predict.py",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "shared",
    "--purpose",
    (
        "Export frozen CYPMol (CjmTH/CYPMol 0c6657d, inhibitor sub0) 512-d molecule-by-isoform "
        "vectors and inhibitor probabilities for all 6896 compounds, mean over three conformers"
    ),
    "--config",
    str(root / "configs/shared/cypmol-embed-v1.json"),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", str(root / ".pixi/envs/cypmol/bin/python"), "scripts/cypmol_embed.py"],
    cwd=root,
    check=True,
)
