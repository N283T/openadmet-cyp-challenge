"""Launch reflected classification leaderboard result decomposition."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
paths = [
    root
    / "artifacts/candidates/classification/20260910T133115130724Z-48208365/predictions.parquet",
    root
    / "artifacts/candidates/classification/20260912T122712113275Z-c0d9d884/predictions.parquet",
    root
    / "artifacts/candidates/classification/20260913T113803083052Z-1322afeb/predictions.parquet",
    root
    / (
        "artifacts/api/classification/20260910T135018498829Z-6408f824/downloads/"
        "dd1e6d0c3ce8a937e47cfc00277d320c7daa864891bf2ba4b3a172a95c25946b/"
        "classification_leaderboard_7rnut467.csv"
    ),
    root
    / (
        "artifacts/api/classification/20260910T135018507465Z-d8cffd06/downloads/"
        "e436029953c626a50c08a8c7200b073c1b85c1eeefea929b015a9b43968a1e24/"
        "classification_leaderboard_k1avqsbk.csv"
    ),
    root
    / (
        "artifacts/api/classification/20260913T113818723423Z-c1045409/downloads/"
        "4ee1b740e6ad43f639f8b2182545c916542bf367efd383ae91a5156d9c8fe137/"
        "classification_leaderboard_2h27jmvh.csv"
    ),
    root
    / (
        "artifacts/api/classification/20260913T113822348087Z-eb1ff1bf/downloads/"
        "5ae8dca377feac816abbdb0ed5dd18a1c2d7f9962f07ad925a394f0ec1e03712/"
        "classification_leaderboard_xyum5p33.csv"
    ),
    root
    / (
        "artifacts/api/classification/20260913T115952558598Z-27dbd92d/downloads/"
        "4a691a64d0ec9a57039386d02ccd40e4a9d40f6765b223a20e56c735e2825e28/"
        "classification_leaderboard_huoq0wk6.csv"
    ),
    root
    / (
        "artifacts/api/classification/20260913T115956532758Z-d04c1f66/downloads/"
        "d0a7146fdcd5153018618894d94208e1fcd7c669d5f75655b6a78095e12d4d3e/"
        "classification_leaderboard_ocp6e_6t.csv"
    ),
    root
    / (
        "artifacts/api/classification/20260913T113815093261Z-d63390f5/downloads/"
        "7670e5e12a7875b5366e5955d11f067d78ac39d98c332ab214d86a778d53b725/"
        "classification_leaderboard_lzvrqs92.csv"
    ),
    root
    / (
        "artifacts/api/classification/20260913T115948535264Z-e392c2cb/downloads/"
        "93c9c376a8329d1e7b87c7dcdfadb7746a4a1868e51936c4e6bfdaceca5c2504/"
        "classification_leaderboard_ojhs6zla.csv"
    ),
]
assert all(path.is_file() for path in paths)
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Decompose the reflected D6-diversity and A4-reranking leaderboard improvement",
    "--config",
    str(root / "configs/classification-round/d6-diverse-lb-analysis-v2.json"),
]
for path in paths:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_lb_result_analysis.py"],
    cwd=root,
    check=True,
)
