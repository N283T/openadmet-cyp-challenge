"""Launch LB-informed fixed-rate rank and blend screening."""

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
    "ensemble_oof": add(
        "artifacts/runs/classification/20260912T121146867140Z-b0a05986/outputs/"
        "ensemble-oof.parquet"
    ),
    "ensemble_test": add(
        "artifacts/runs/classification/20260912T121146867140Z-b0a05986/outputs/"
        "ensemble-test.parquet"
    ),
    "caruana_oof": add(
        "artifacts/runs/classification/20260912T141520708830Z-97df4c71/outputs/oof.parquet"
    ),
    "caruana_test": add(
        "artifacts/runs/classification/20260912T141520708830Z-97df4c71/outputs/"
        "deployment.parquet"
    ),
    "weighted_oof": add(
        "artifacts/runs/classification/20260913T102143243794Z-b9ebbcba/outputs/oof.parquet"
    ),
    "weighted_test": add(
        "artifacts/runs/classification/20260913T102143243794Z-b9ebbcba/outputs/"
        "deployment.parquet"
    ),
    "old_submission": add(
        "artifacts/candidates/classification/20260910T133115130724Z-48208365/predictions.parquet"
    ),
    "current_submission": add(
        "artifacts/candidates/classification/20260912T122712113275Z-c0d9d884/predictions.parquet"
    ),
    "aliases": add(
        "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/aliases.parquet"
    ),
    "lb_summary": add(
        "artifacts/runs/classification/20260913T100729603583Z-d7f1610c/outputs/summary.json"
    ),
    "our_lb": add(
        "artifacts/runs/classification/20260913T100729603583Z-d7f1610c/outputs/"
        "our-lb-operating-points.csv"
    ),
    "leaderboard_points": add(
        "artifacts/runs/classification/20260913T100729603583Z-d7f1610c/outputs/"
        "leaderboard-operating-points.csv"
    ),
    "endpoints": ["CYP2D6_is_TDI", "CYP3A4_is_TDI"],
}
config = root / "configs/classification-round/lb-rank-screen-v1.json"
config.write_text(json.dumps(cfg, indent=2) + "\n")
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Screen fixed LB-derived positive rates and rank blends on saved TDI predictions",
    "--config",
    str(config),
]
for path in paths:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_lb_rank_screen.py"],
    cwd=root,
    check=True,
)
