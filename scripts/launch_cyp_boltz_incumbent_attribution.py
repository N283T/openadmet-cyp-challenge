"""Add a no-trunk reselection control without refitting any saved readouts."""

import json
import subprocess
import sys
from pathlib import Path

from launch_cyp_boltz_incumbent_suite import recorded


def main():
    root = Path(__file__).resolve().parents[1]
    directory = Path(sys.argv[1]).resolve()
    previous_id = (directory / "regression-placement-run.txt").read_text().strip()
    previous = root / "artifacts/runs/regression" / previous_id
    manifest = json.loads((previous / "manifest.json").read_text())
    cfg = json.loads((previous / "config.json").read_text())
    cfg["attribution_control"] = "Same bank reselection algorithm without trunk"
    inputs = {Path(item["snapshot"]).name: Path(item["snapshot"]) for item in manifest["inputs"]}
    identifier = recorded(
        directory, f"regression-incumbent-attribution-{previous_id.rsplit('-', 1)[-1]}",
        "regression", cfg, inputs,
        "evaluate_cyp_boltz_incumbent.py",
    )
    (directory / "previous-regression-placement-run.txt").write_text(previous_id + "\n")
    (directory / "regression-placement-run.txt").write_text(identifier + "\n")
    subprocess.run(
        ["pixi", "run", "python", "scripts/report_cyp_boltz_incumbent.py", str(directory)],
        cwd=root, check=True,
    )


if __name__ == "__main__":
    main()
