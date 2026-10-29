"""Launch the pre-registration of the next D6 call-count probe."""

import argparse
import subprocess
from pathlib import Path

DATA = "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
SEP15_BLIND = "artifacts/runs/classification/20260924T045344348206Z-7e2db311/inputs"
TESTLIKE = "artifacts/tracked/cyp-d6-testlike-pool-20261002T215345688167Z-7920b49a"
PAYLOAD = "artifacts/submissions/classification/20261002T223229553466Z-a08301d2"
INPUTS = {
    "payload-oct03.parquet": f"{PAYLOAD}/predictions.parquet",
    "sep15-blind-scores.parquet": f"{SEP15_BLIND}/blind-scores.parquet",
    "pool-blind-scores.parquet": f"{TESTLIKE}/test_scores.parquet",
    "test_mapping.parquet": f"{DATA}/test_mapping.parquet",
    "testlike-pool-oof.parquet": f"{TESTLIKE}/oof_scores.parquet",
}


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config", default="configs/classification/d6-probe-preregistration-v1.json"
    )
    args = parser.parse_args()
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        "Pre-register the next D6 call-count probe: reproduce the deployed L4 ordering, "
        "count changed calls for each candidate count against the October 3 payload, and "
        "write the public-metric prediction as a function of the marginal precision of the "
        "added ranks; nothing fitted, no candidate, no prepared submission, no upload",
        "--config",
        str(root / args.config),
    ]
    for alias, relative in INPUTS.items():
        path = root / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        command.extend(["--input", f"{alias}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/evaluate_cyp_d6_probe_preregistration.py"],
        cwd=root,
        check=True,
    )


if __name__ == "__main__":
    main()
