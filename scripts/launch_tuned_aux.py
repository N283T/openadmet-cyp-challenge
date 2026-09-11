"""Refit a preselected HPO recipe for the untouched outer-fold-zero readout."""

import argparse
import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("--search-run", required=True)
args = parser.parse_args()
shared = root / "artifacts/runs/shared"
selected_path = shared / args.search_run / "outputs/selected.json"
selected = json.loads(selected_path.read_text())
for other in range(1, 5):
    cfg = {
        **selected["config"],
        "tuning_only": False,
        "excluded_folds": [0, other],
        "outer_fold": 0,
        "selection_run": args.search_run,
        "selection_trial": selected["trial"],
    }
    p = root / f"configs/parent-cv/tuned-aux-0-{other}.json"
    with p.open("x") as stream:
        json.dump(cfg, stream, indent=2)
    cmd = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "shared",
        "--purpose",
        f"Evaluate frozen Rustuna choice in nested auxiliary pair 0,{other}",
        "--config",
        str(p),
    ]
    for path in [
        shared / "20260910T100337258849Z-927eac1c/outputs/targets.parquet",
        shared / "20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
        root
        / "artifacts/runs/regression/20260911T100134004950Z-c7a3f2bf/outputs/assignments.parquet",
        shared / "20260910T012612848389Z-bfeed434/outputs/chemeleon_mp.pt",
        shared / "20260910T012612848389Z-bfeed434/outputs/provenance.json",
        selected_path,
    ]:
        cmd += ["--input", str(path)]
    subprocess.run(
        cmd + ["--", "python", "-m", "cyp_challenge.log2fc_pretrain"], check=True, cwd=root
    )
