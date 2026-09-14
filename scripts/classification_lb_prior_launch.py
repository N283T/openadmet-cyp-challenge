"""Launch the leaderboard-prior and threshold analysis with frozen inputs."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
paths = []


def add(path):
    candidate = root / path
    assert candidate.is_file(), candidate
    paths.append(candidate)
    return f"{len(paths) - 1:03d}-{candidate.name}"


cfg = {
    "old_submission": add(
        "artifacts/candidates/classification/20260910T133115130724Z-48208365/predictions.parquet"
    ),
    "current_submission": add(
        "artifacts/candidates/classification/20260912T122712113275Z-c0d9d884/predictions.parquet"
    ),
    "old_d6_leaderboard": add(
        "artifacts/api/classification/20260910T135018498829Z-6408f824/downloads/"
        "dd1e6d0c3ce8a937e47cfc00277d320c7daa864891bf2ba4b3a172a95c25946b/"
        "classification_leaderboard_7rnut467.csv"
    ),
    "old_a4_leaderboard": add(
        "artifacts/api/classification/20260910T135018507465Z-d8cffd06/downloads/"
        "e436029953c626a50c08a8c7200b073c1b85c1eeefea929b015a9b43968a1e24/"
        "classification_leaderboard_k1avqsbk.csv"
    ),
    "current_d6_leaderboard": add(
        "artifacts/api/classification/20260912T124752635013Z-591d549f/downloads/"
        "fe590b0b6a0b17b0ddcfc5c2831386ac5f5fb061cc44e64bb5f8124ac77a75f7/"
        "classification_leaderboard_wairi7uv.csv"
    ),
    "current_a4_leaderboard": add(
        "artifacts/api/classification/20260912T124752634908Z-3ddcd27d/downloads/"
        "682a3fcf50b45899acbe8f4fe372db064971407e689d780b1f9c92fb8bbedd14/"
        "classification_leaderboard_elbwovr6.csv"
    ),
    "ensemble_oof": add(
        "artifacts/runs/classification/20260912T121146867140Z-b0a05986/outputs/"
        "ensemble-oof.parquet"
    ),
    "ensemble_inner": add(
        "artifacts/runs/classification/20260912T121146867140Z-b0a05986/outputs/"
        "ensemble-inner.parquet"
    ),
    "ensemble_test": add(
        "artifacts/runs/classification/20260912T121146867140Z-b0a05986/outputs/"
        "ensemble-test.parquet"
    ),
    "aliases": add(
        "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/aliases.parquet"
    ),
    "official_evaluation": add(
        "artifacts/runs/classification/20260912T130054345728Z-e960da8e/inputs/"
        "012-evaluate_predictions.py"
    ),
    "endpoints": ["CYP2D6_is_TDI", "CYP3A4_is_TDI"],
    "current_methods": {
        "CYP2D6_is_TDI": "all_plain_pfn_graph_mean",
        "CYP3A4_is_TDI": "lgbm_balanced_mean",
    },
    "all_negative_username": "cyp450_reg",
    "our_username": "N283T",
}
config = root / "configs/classification-round/lb-prior-v1.json"
config.write_text(json.dumps(cfg, indent=2) + "\n")

command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Infer hidden scored-class priors from public LB and test prior-aware TDI thresholds",
    "--config",
    str(config),
]
for path in paths:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_lb_prior_analysis.py"],
    cwd=root,
    check=True,
)
