"""Launch the statistical post-mortem of the 2026-10-03 D6 Live MCC gain."""

import argparse
import subprocess
from pathlib import Path

API = "artifacts/api/classification"
DATA = "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
POOL_RUN = "artifacts/runs/classification/20261002T215345688167Z-7920b49a/outputs"
SEP15_BLIND = "artifacts/runs/classification/20260924T045344348206Z-7e2db311/inputs"
CONCORDANCE = "artifacts/runs/classification/20261002T221918697661Z-2d8d20f5/outputs"
DIAGNOSTIC = "artifacts/runs/classification/20261002T223515148996Z-ccb05cf3/outputs"
SUBMISSIONS = "artifacts/submissions/classification"

# Saved D6 Live leaderboard snapshots, one per distinct N283T submission timestamp.
LEADERBOARDS = {
    "lb_sep13.csv": "20260913T115952558598Z-27dbd92d",
    "lb_sep14.csv": "20260914T021850504300Z-f8bc9ee4",
    "lb_sep15.csv": "20260915T130816907483Z-fa8e6053",
    "lb_sep18.csv": "20260918T101216302522Z-4a8b944c",
    "lb_sep23.csv": "20260923T183846404373Z-b2478be3",
    "lb_oct01.csv": "20261001T142809834162Z-478c049d",
    "lb_oct02.csv": "20261002T125228445980Z-10056d8d",
    "lb_oct03.csv": "20261003T004710939132Z-27ab5110",
}
# The MA and A4 boards of the same 2026-10-03 download batch, for provenance.
CONTEXT = {
    "lb_oct03_ma.csv": "20261003T004649176162Z-9e8eb70d",
    "lb_oct03_a4.csv": "20261003T004717577185Z-228dc8dc",
}
PAYLOADS = {
    "payload_sep15.parquet": "20260915T130004513006Z-75b25070",
    "payload_oct03.parquet": "20261002T223229553466Z-a08301d2",
}
PURPOSE = (
    "D6 Live gain post-mortem: reconstruct the scored subset and every payload confusion "
    "matrix from the published bootstrap aggregates, size the Sep 15 to Oct 3 change as row "
    "outcomes, resample an unskilled 17-call swap null, attribute the new vector between the "
    "deployed score and pool_lgbm, and extend the concordance table by the eighth "
    "transition; read-only, no candidate or upload"
)


def leaderboard_csv(root, snapshot):
    matches = sorted((root / API / snapshot / "downloads").glob("*/*.csv"))
    if len(matches) != 1:
        raise FileNotFoundError(f"{snapshot}: expected one leaderboard CSV, found {len(matches)}")
    return matches[0]


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/classification/d6-live-gain-v1.json")
    args = parser.parse_args()

    inputs = {alias: leaderboard_csv(root, s) for alias, s in LEADERBOARDS.items()}
    inputs.update({alias: leaderboard_csv(root, s) for alias, s in CONTEXT.items()})
    inputs.update(
        {
            alias: root / f"{SUBMISSIONS}/{run}/predictions.parquet"
            for alias, run in PAYLOADS.items()
        }
    )
    inputs.update(
        {
            "test_scores.parquet": root / POOL_RUN / "test_scores.parquet",
            "sep15_blind_scores.parquet": root / SEP15_BLIND / "blind-scores.parquet",
            "test_mapping.parquet": root / DATA / "test_mapping.parquet",
            "call_diagnostics.parquet": root / DIAGNOSTIC / "call_diagnostics.parquet",
            "changed_calls.csv": root / DIAGNOSTIC / "changed_calls.csv",
            "transitions.csv": root / CONCORDANCE / "transitions.csv",
        }
    )
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        PURPOSE,
        "--config",
        str(root / args.config),
        "--input-mode",
        "hardlink",
    ]
    for alias, path in inputs.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        command.extend(["--input", f"{alias}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/evaluate_cyp_d6_live_gain.py"], cwd=root, check=True
    )


if __name__ == "__main__":
    main()
