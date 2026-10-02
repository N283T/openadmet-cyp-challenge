"""Record exact target-column composition for local inspection only."""

import subprocess
from pathlib import Path

from cyp_challenge.database import connect
from cyp_challenge.submission import last_sent


def main():
    root = Path(__file__).resolve().parents[1]
    submissions = root / "artifacts/submissions/classification"
    expected = submissions / "20260926T150948735129Z-6fc9d05b/predictions.parquet"
    with connect() as conn:
        _, current_last = last_sent(conn, "classification")
    if current_last != expected:
        raise ValueError("Last sent payload changed; update the review configuration")
    inputs = {
        "d6-probe.parquet": submissions / "20261001T135846483051Z-64386ab4/predictions.parquet",
        "best-control.parquet": submissions / "20260924T050705396203Z-10e1d3ad/predictions.parquet",
        "last-sent.parquet": current_last,
        "reference.csv": root / "data/raw/cyp-challenge-TEST-BLINDED.csv",
    }
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        "Prepare exact reviewed D6 AP probe with best-scored A4 restored; no upload",
        "--config",
        str(root / "configs/boltz-trunk/d6-ap-probe-a4-best-restoration-v1.json"),
    ]
    for alias, path in inputs.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        command.extend(["--input", f"{alias}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/package_cyp_boltz_d6_restore_a4.py"],
        cwd=root,
        check=True,
    )


if __name__ == "__main__":
    main()
