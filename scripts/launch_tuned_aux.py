"""Refit a preselected HPO recipe for the untouched outer-fold-zero readout."""

import argparse
import itertools
import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument("--search-run", required=True)
parser.add_argument(
    "--remaining", action="store_true", help="Fit the six pairs not containing fold zero"
)
args = parser.parse_args()
shared = root / "artifacts/runs/shared"
selected_path = shared / args.search_run / "outputs/selected.json"
selected = json.loads(selected_path.read_text())
for first, other in (
    itertools.combinations(range(1, 5), 2) if args.remaining else [(0, i) for i in range(1, 5)]
):
    cfg = {
        **selected["config"],
        "tuning_only": False,
        "excluded_folds": [first, other],
        "outer_fold": first,
        "selection_run": args.search_run,
        "selection_trial": selected["trial"],
    }
    p = root / f"configs/parent-cv/tuned-aux-{first}-{other}.json"
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
        f"Fit adopted Rustuna choice in nested auxiliary pair {first},{other}",
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
