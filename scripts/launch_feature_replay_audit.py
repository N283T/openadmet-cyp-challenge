"""Record replay inputs without reading mutable assets inside the audit."""

import json
import subprocess
from pathlib import Path

base = Path("artifacts/runs")
files = dict(
    local=Path("data/raw/cyp-challenge-TRAIN_Emax.csv"),
    remote=base / "regression/20260911T141511080701Z-fd93c411/outputs/cyp-challenge-TRAIN_Emax.csv",
    compounds=base / "shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    standardization=base / "shared/20260909T132734119233Z-a7063ccd/outputs/standardization.parquet",
    embeddings=base / "shared/20260910T012630927160Z-e9576045/outputs/embeddings.parquet",
    checkpoint=base / "shared/20260910T012612848389Z-bfeed434/outputs/chemeleon_mp.pt",
    metadata=base / "shared/20260911T123506505890Z-15c68a1d/outputs/metadata.json",
    best=base / "shared/20260911T123506505890Z-15c68a1d/outputs/best.pt",
)
inputs = list(files.values())
config = {k: f"{i:03d}-{v.name}" for i, (k, v) in enumerate(files.items())}
config["predictions"] = []
for f in sorted(
    (base / "shared/20260911T123506505890Z-15c68a1d/outputs").glob("predictions-*.parquet")
):
    config["predictions"].append(f"{len(inputs):03d}-{f.name}")
    inputs.append(f)
p = Path("configs/regression/feature-replay-audit.json")
p.write_text(json.dumps(config, indent=2) + "\n")
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Replay frozen embeddings and tuned auxiliary predictions from SMILES; inspect Emax update",
    "--config",
    str(p),
]
for f in inputs:
    cmd += ["--input", str(f)]
subprocess.run(cmd + ["--", "python", "scripts/audit_feature_replay.py"], check=True)
