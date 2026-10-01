"""Record authorized two-target regression adoption and freeze a candidate."""

import json
import subprocess
from datetime import UTC, datetime
from pathlib import Path

from launch_cyp_boltz_incumbent_suite import recorded


def main():
    root = Path(__file__).resolve().parents[1]
    source_id = "20261001T111723965954Z-04ff278c"
    runs = root / "artifacts/runs/regression"
    source = runs / source_id
    cfg = json.loads((source / "config.json").read_text())
    cfg["targets"] = {t: cfg["targets"][t] for t in ["cyp1a2", "cyp2d6"]}
    cfg["best_submission"] = "20260926T090548323924Z-46188564"
    cfg["last_submission"] = "20260926T140601457029Z-bd3a1e40"
    cfg["multi_analysis"] = "20261001T113601716233Z-bf06a273"
    cfg["deployment_policy"] = (
        "1A2 multiple replacement; D6 single trunk bank extension; final best moments"
    )
    inputs = {
        Path(item["snapshot"]).name: Path(item["snapshot"])
        for item in json.loads((source / "manifest.json").read_text())["inputs"]
        if Path(item["snapshot"]).name.startswith(("cyp1a2-", "cyp2d6-", "d6-base"))
    }
    inputs.update(
        {
            "best.csv": root
            / "artifacts/submissions/regression"
            / cfg["best_submission"]
            / "predictions.csv",
            "last.csv": root
            / "artifacts/submissions/regression"
            / cfg["last_submission"]
            / "predictions.csv",
            "reference.csv": root / "data/raw/cyp-challenge-TEST-BLINDED.csv",
            "aliases.parquet": root
            / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/aliases.parquet",
            "atomic-test.parquet": runs
            / "20260918T035106389077Z-b97ec048/outputs/atomic-test-predictions.parquet",
            "d6-core-test.parquet": runs
            / "20260912T163036856426Z-bf2e35dc/outputs/test-predictions.parquet",
            "d6-robust-test.parquet": runs
            / "20260926T090408836539Z-239369e6/outputs/d6-test-comparison.parquet",
        }
    )
    directory = (
        root
        / "artifacts/boltz-trunk/experiments"
        / (datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ") + "-adoption")
    )
    directory.mkdir()
    identifier = recorded(
        directory,
        "package-trunk-regression",
        "regression",
        cfg,
        inputs,
        "package_cyp_boltz_regression.py",
    )
    print(f"package_run={identifier}", flush=True)
    subprocess.run(
        [
            "pixi",
            "run",
            "cyp",
            "candidate",
            "--track",
            "regression",
            "--file",
            str(runs / identifier / "outputs/predictions.csv"),
            "--reference",
            str(inputs["reference.csv"]),
            "--run",
            identifier,
            "--notes",
            "User-adopted Boltz 1A2 multiple replacement and D6 single-trunk extension; "
            "C9/A4 exact best submitted vectors; preserve best calibration moments.",
        ],
        cwd=root,
        check=True,
    )


if __name__ == "__main__":
    main()
