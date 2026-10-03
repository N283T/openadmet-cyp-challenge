"""Launch the read-only changed-call diagnostic of the two prepared D6 gate candidates."""

import argparse
import json
import subprocess
from pathlib import Path

RUNS = "artifacts/runs"
DATA = f"{RUNS}/shared/20260909T132734119233Z-a7063ccd/outputs"
SUBMISSIONS = "artifacts/submissions/classification"
SEP15_SUBMISSION = "20260924T050705396203Z-10e1d3ad"
PURPOSE = (
    "Read-only diagnostic of the prepared D6 gate candidates L4 (a08301d2) and L2 "
    "(16026b6c): legacy, pool_lgbm and combined ranks, Morgan-gate p_active and "
    "similarity-weighted top-10 labelled D6 neighbour positive rate for retained, "
    "removed and added calls versus the Sep 15 payload; no candidate or upload"
)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config", default="configs/classification/d6-gate-candidate-diagnostic-v1.json"
    )
    args = parser.parse_args()
    cfg = json.loads((root / args.config).read_text())
    inputs = {
        "sep15_payload.parquet": f"{SUBMISSIONS}/{SEP15_SUBMISSION}/predictions.parquet",
        "test_scores.parquet": f"{RUNS}/classification/{cfg['gate_run']}"
        "/outputs/test_scores.parquet",
        "compounds.parquet": f"{DATA}/compounds.parquet",
        "measurements.parquet": f"{DATA}/measurements.parquet",
    }
    for spec in cfg["candidates"].values():
        inputs[spec["scores_input"]] = (
            f"{RUNS}/classification/{spec['package_run']}/outputs/test_scores.parquet"
        )
        inputs[spec["prepared_input"]] = f"{SUBMISSIONS}/{spec['prepared_id']}/predictions.parquet"
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
    for alias, relative in inputs.items():
        path = root / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        command.extend(["--input", f"{alias}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/diagnose_cyp_d6_gate_candidates.py"],
        cwd=root,
        check=True,
    )


if __name__ == "__main__":
    main()
