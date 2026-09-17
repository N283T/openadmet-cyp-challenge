"""Launch recorded curation of pinned OpenADMET metabolism sources."""

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FETCH = ROOT / "artifacts/runs/shared/20260916T072328099350Z-fa03aef6/outputs"
BASE = ROOT / "artifacts/runs/shared/20260912T055809060996Z-06088efb/outputs"
CHALLENGE_CONNECTIVITY = (
    ROOT
    / "artifacts/runs/shared/20260911T052745709989Z-8792f422/outputs/challenge-connectivity.parquet"
)
CONFIG = ROOT / "configs/external-source/openadmet-metabolism-curation-v1.json"

inputs = {
    name: FETCH / name
    for name in [
        "octant-inhibition.tsv",
        "octant-reactivity.tsv",
        "octant-mass-spec.tsv",
        "expansionrx-raw.csv",
        "pxr-train.csv",
        "pxr-phase1.csv",
        "pxr-phase2.csv",
    ]
}
inputs.update(
    {
        "external-compounds.parquet": BASE / "external-compounds.parquet",
        "observations.parquet": BASE / "observations.parquet",
        "challenge-connectivity.parquet": CHALLENGE_CONNECTIVITY,
    }
)
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "shared",
    "--purpose",
    "Curate pinned OpenADMET Octant, ExpansionRx, and PXR metabolism data; "
    "deduplicate older Biogen structures and exclude all challenge connectivity",
    "--config",
    str(CONFIG),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/openadmet_metabolism_curate.py"],
    cwd=ROOT,
    check=True,
)
