"""Refit the original CheMeleon auxiliary recipe on the new series groups."""

import itertools
import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
shared = root / "artifacts/runs/shared"
split = "20260911T100134004950Z-c7a3f2bf"
for a, b in itertools.combinations(range(5), 2):
    cfg = json.loads((root / f"configs/stacking/aux-{a}-{b}.json").read_text())
    cfg.update(
        split_run=split,
        outer_strategy="series",
        purpose="Parent CV auxiliary: exclude outer/recipient groups including parents",
    )
    path = root / f"configs/parent-cv/aux-{a}-{b}.json"
    path.parent.mkdir(exist_ok=True)
    with path.open("x") as f:
        json.dump(cfg, f, indent=2)
    paths = [
        shared / "20260910T100337258849Z-927eac1c/outputs/targets.parquet",
        shared / "20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
        root / f"artifacts/runs/regression/{split}/outputs/assignments.parquet",
        shared / "20260910T012612848389Z-bfeed434/outputs/chemeleon_mp.pt",
        shared / "20260910T012612848389Z-bfeed434/outputs/provenance.json",
    ]
    cmd = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "shared",
        "--purpose",
        f"Refit CheMeleon log2fc auxiliary excluding series groups in folds {a},{b}",
        "--config",
        str(path),
    ]
    for p in paths:
        cmd += ["--input", str(p)]
    subprocess.run(
        cmd + ["--", "python", "-m", "cyp_challenge.log2fc_pretrain"], check=True, cwd=root
    )
