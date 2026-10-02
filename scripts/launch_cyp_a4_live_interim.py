"""Launch the A4 Live-versus-Interim leaderboard separation analysis."""

import argparse
import subprocess
from pathlib import Path

API = "artifacts/api/classification"
INPUTS = {
    "live-a4.csv": f"{API}/20260926T001919860383Z-58e63f72/downloads/"
    "c8c241e5e0f8d75368ec9dc64ec496b2522ca11f42bce57c5cdc27eec9b1af44/"
    "classification_leaderboard_dvvui32a.csv",
    "live-d6.csv": f"{API}/20260926T001915310713Z-c9c2c8e5/downloads/"
    "a8b8cae2ee05a612a52a006b72eca3bf2457215b701ba5b43844eb42df51327a/"
    "classification_leaderboard_6xhyv75f.csv",
    "live-ma.csv": f"{API}/20260926T001910749563Z-e5ea6039/downloads/"
    "b24d14b78463ad019f6797961e5ae0ad603fd41f8b4f3baa605470c92ff82425/"
    "classification_leaderboard_n9bbauf6.csv",
    "interim-a4.csv": f"{API}/20260926T001951338594Z-19f44815/downloads/"
    "544893776cdc0514d72cfb303b5e87e7968a3334594fa8dd98ebd67b1de9c35d/"
    "classification_leaderboard_lagiiqu7.csv",
    "interim-d6.csv": f"{API}/20260926T001946955634Z-044a1693/downloads/"
    "91e6ea1225965089fa873fc0ee5322519a2c9fe9f75556f7298da75eb34618aa/"
    "classification_leaderboard_e8upvwm3.csv",
    "interim-ma.csv": f"{API}/20260926T001826505043Z-088d559e/downloads/"
    "65dd8e62ae65a46e37b6e754c30ae1ec21cc0906b4d304bd7ff947a8e5ea1d87/"
    "classification_leaderboard_k64ru3hh.csv",
    "snapshot-manifest.json": "artifacts/tracked/leaderboard-snapshot-20260926/manifest.json",
    "history-a4-20260910.csv": f"{API}/20260910T135018507465Z-d8cffd06/downloads/"
    "e436029953c626a50c08a8c7200b073c1b85c1eeefea929b015a9b43968a1e24/"
    "classification_leaderboard_k1avqsbk.csv",
    "history-a4-20260912.csv": f"{API}/20260912T124752634908Z-3ddcd27d/downloads/"
    "682a3fcf50b45899acbe8f4fe372db064971407e689d780b1f9c92fb8bbedd14/"
    "classification_leaderboard_elbwovr6.csv",
    "history-a4-20260913.csv": f"{API}/20260913T115956532758Z-d04c1f66/downloads/"
    "d0a7146fdcd5153018618894d94208e1fcd7c669d5f75655b6a78095e12d4d3e/"
    "classification_leaderboard_ocp6e_6t.csv",
    "history-a4-20260914.csv": f"{API}/20260914T021854541658Z-77f4b6cd/downloads/"
    "65991b188f12a68012e2e4dd1bec400237a3bde20cfb83863ada96b18ad66513/"
    "classification_leaderboard_6elk3t3j.csv",
    "history-a4-20260915.csv": f"{API}/20260915T130816924118Z-19f008d0/downloads/"
    "4624e63936755de53709bad35c4e47e22666e71693f8d9c23cc404b626c02bed/"
    "classification_leaderboard_n_8cx0gz.csv",
    "history-a4-20260918.csv": f"{API}/20260918T101216302525Z-527127ee/downloads/"
    "d4261f8ba291e9c0c42515c256e7c7920449748049e0004e5df95599d51fc152/"
    "classification_leaderboard_h19mu01y.csv",
    "history-a4-20260923.csv": f"{API}/20260923T181924882069Z-ed5fddd6/downloads/"
    "4857fd79eec032923d8c1985500e88b88400751b54bc42a12aab263e937ae463/"
    "classification_leaderboard_gwfjyw0y.csv",
    "history-a4-20260924.csv": f"{API}/20260924T063722820119Z-e5e2f1a3/downloads/"
    "1f87c2740701cc7562d80e1786f856486cadcdc4a74103ce431c2bfa463d664e/"
    "classification_leaderboard_6i87hvbf.csv",
    "history-a4-20260926.csv": f"{API}/20260926T152036592294Z-8c6a9959/downloads/"
    "c010d39bd7d43f1e3805301f1a7985e84204ad1abe67a580ab90925b3566e736/"
    "classification_leaderboard_l91olpq7.csv",
    "history-a4-20261001.csv": f"{API}/20261001T142809834162Z-9ada1e4e/downloads/"
    "3c1f342e1d9b794757920c1ac69eeef96233bb9619b2189ea7d75b2db9e151ec/"
    "classification_leaderboard_83icvy5h.csv",
    "history-a4-20261002.csv": f"{API}/20261002T125228445972Z-6fdebba0/downloads/"
    "8ba820f88522d216d794eb2989ab56f003bce47727714c83ae3ef4c6bfb446ca/"
    "classification_leaderboard_x06bdmv5.csv",
}


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/classification/a4-live-interim-v1.json")
    args = parser.parse_args()
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        "A4 TDI: compare the September 26 Live and Interim CYP3A4 boards, express the "
        "Live rank-1 lead and the Interim deficit in bootstrap std units, infer scored "
        "positive counts and call fractions, and list the N283T A4 Live history",
        "--config",
        str(root / args.config),
        "--input-mode",
        "hardlink",
    ]
    for alias, relative in INPUTS.items():
        path = root / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        command.extend(["--input", f"{alias}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/evaluate_cyp_a4_live_interim.py"],
        cwd=root,
        check=True,
    )


if __name__ == "__main__":
    main()
