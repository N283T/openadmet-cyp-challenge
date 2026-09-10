"""Launch OOF audit and register completed stacked readouts (no model fitting)."""

import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    log = ROOT / "artifacts/launches/stacking-readout-20260910.log"
    runs = re.findall(r"run_id=(\S+)", log.read_text())
    assert len(runs) == len(set(runs)) == 5
    cmd = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "shared",
        "--purpose",
        "Audit PCA256 plus nested predicted log2fc TabPFN OOF",
        "--config",
        "configs/readout/audit.json",
        "--input",
        "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
        "--input",
        "artifacts/metric-reference/custom_scoring_functions.py",
    ]
    for run in runs:
        for name in ["oof.parquet", "summary.json"]:
            cmd += ["--input", str(ROOT / "artifacts/runs/shared" / run / "outputs" / name)]
    subprocess.run(cmd + ["--", "python", "scripts/audit_readouts.py"], cwd=ROOT, check=True)
    cmd = ["pixi", "run", "python", "scripts/register_stacking.py"]
    for run in runs:
        cmd += ["--run", run]
    subprocess.run(cmd, cwd=ROOT, check=True)


if __name__ == "__main__":
    main()
