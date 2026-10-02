"""Record a blind-call review using retained scores; no remote submission route."""

import argparse
import json
import subprocess
from pathlib import Path

from cyp_challenge.database import connect
from cyp_challenge.submission import last_sent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("selection_run")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    runs = root / "artifacts/runs/classification"
    with connect() as conn:
        sent_at, last = last_sent(conn, "classification")
    if last is None:
        raise ValueError("No recorded last sent classification payload")
    last_id = last.parent.name
    inputs = {
        "decision.json": runs / args.selection_run / "outputs/decision.json",
        "control-scores.parquet": runs
        / "20260915T125939829893Z-1d921ffc/outputs/test-scores.parquet",
        "trunk-predictions.parquet": runs
        / "20261001T023655009132Z-91096754/outputs/predictions.parquet",
        "last-sent.parquet": last,
        "performance-control.parquet": root
        / "artifacts/submissions/classification/20260924T050705396203Z-10e1d3ad"
        / "predictions.parquet",
        "aliases.parquet": root
        / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/aliases.parquet",
        "reference.csv": root / "data/raw/cyp-challenge-TEST-BLINDED.csv",
        "neighbor-diagnostics.parquet": runs
        / "20260923T190811462293Z-4f4a92f2/outputs/threeway-test.parquet",
    }
    cfg = {
        "selection_run": args.selection_run,
        "last_sent_submission": last_id,
        "last_sent_at": sent_at.isoformat(),
        "performance_control_submission": "20260924T050705396203Z-10e1d3ad",
        "scope": (
            "User authorized inspection only: D6-only blind review, "
            "A4 last-sent column preserved; no upload"
        ),
    }
    config = root / "configs/boltz-trunk" / f"d6-ap-probe-package-{args.selection_run}.json"
    if config.exists():
        if json.loads(config.read_text()) != cfg:
            raise ValueError("Existing config differs")
    else:
        with config.open("x") as stream:
            json.dump(cfg, stream, indent=2)
            stream.write("\n")
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        "Inspect 750 blind D6 AP-trunk rank-blend calls; freeze only, no upload",
        "--config",
        str(config),
    ]
    for alias, path in inputs.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        command.extend(["--input", f"{alias}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/package_cyp_boltz_d6_ap_probe.py"], cwd=root, check=True
    )


if __name__ == "__main__":
    main()
