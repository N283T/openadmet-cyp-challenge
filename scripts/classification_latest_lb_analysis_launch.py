"""Launch the recorded latest classification leaderboard analysis."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
api = root / "artifacts/api/classification"
submissions = root / "artifacts/submissions/classification"
runs = root / "artifacts/runs/classification"
inputs = [
    api / "20260915T130013900061Z-5aa75bfd/downloads"
    / "415b2407bcf1ab0b8cfde1eee30b0a322ecdbb563fb8ac42cefe272fc55ba7f0"
    / "classification_leaderboard_u5msmx90.csv",
    api / "20260915T130017486614Z-f261b37c/downloads"
    / "98f9635bccce98ff799c3b891c3a8007dcf3b5ecb9c6f361b7e05980dd5d20d9"
    / "classification_leaderboard_n71zlaf1.csv",
    api / "20260915T130021049578Z-4cb4e828/downloads"
    / "62621a18f75222091ea02e70964842b7e30a47e594237ec0baacbb22d61ae310"
    / "classification_leaderboard_tkhrqrpx.csv",
    api / "20260915T130749950536Z-4015ea99/downloads"
    / "0f28c67f83bc39645d17146deff7150a63de2e0e4fc39b3bdbf978b412ddfebb"
    / "classification_leaderboard_8782ot6k.csv",
    api / "20260915T130816907483Z-fa8e6053/downloads"
    / "4938e4a36b83370a323761258778fddc2270a93d054cf6a955b12dd789e63e38"
    / "classification_leaderboard_uorsc9xo.csv",
    api / "20260915T130816924118Z-19f008d0/downloads"
    / "4624e63936755de53709bad35c4e47e22666e71693f8d9c23cc404b626c02bed"
    / "classification_leaderboard_n_8cx0gz.csv",
    submissions / "20260914T020759952207Z-37dd7bd8/predictions.parquet",
    submissions / "20260915T130004513006Z-75b25070/predictions.parquet",
    runs / "20260915T125001009276Z-d5adff22/outputs/bootstrap.csv",
    runs / "20260915T113718060068Z-b3a7120e/outputs/metrics.csv",
]
for path in inputs:
    if not path.is_file():
        raise FileNotFoundError(path)
command = [
    "pixi", "run", "cyp", "run", "--track", "classification",
    "--purpose", "Analyze reflected D6 mean-block plus A4 equal-seven leaderboard result",
    "--config", str(root / "configs/analysis/classification-latest-lb-v1.json"),
]
for path in inputs:
    command += ["--input", str(path)]
subprocess.run(
    command + ["--", "python", "scripts/classification_latest_lb_analysis.py"],
    cwd=root,
    check=True,
)
