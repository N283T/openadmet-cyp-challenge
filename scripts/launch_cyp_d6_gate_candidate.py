"""Launch the packaging run of one legacy x gate D6 candidate variant.

The base config names a default `variant`; `--variant` writes a resolved copy with
that field replaced, which the run wrapper snapshots as the run config.
"""

import argparse
import json
import subprocess
import tempfile
from pathlib import Path

from cyp_challenge.database import connect
from cyp_challenge.submission import last_sent

RUNS = "artifacts/runs"
SUBMISSIONS = "artifacts/submissions/classification"


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/classification/d6-gate-candidate-v1.json")
    parser.add_argument("--variant")
    args = parser.parse_args()
    cfg = json.loads((root / args.config).read_text())
    if args.variant:
        cfg["variant"] = args.variant
    if cfg["variant"] not in cfg["variants"]:
        raise ValueError(f"Unknown variant: {cfg['variant']}")
    with connect() as conn:
        _, last = last_sent(conn, "classification")
    if last is None or last.parent.name != cfg["last_sent_submission"]:
        raise ValueError(f"Last sent classification payload is {last}, not the configured one")
    inputs = {
        "sep15_blind_scores.parquet": f"{RUNS}/classification/{cfg['legacy_run']}"
        "/inputs/blind-scores.parquet",
        "test_scores.parquet": f"{RUNS}/classification/{cfg['score_run']}"
        "/outputs/test_scores.parquet",
        "blind_candidates.parquet": f"{RUNS}/classification/{cfg['selection_run']}"
        "/outputs/blind_candidates.parquet",
        "sep15_payload.parquet": f"{SUBMISSIONS}/{cfg['sep15_submission']}/predictions.parquet",
        "last_sent.parquet": f"{SUBMISSIONS}/{cfg['last_sent_submission']}/predictions.parquet",
        "test_mapping.parquet": f"{RUNS}/shared/20260909T132734119233Z-a7063ccd"
        "/outputs/test_mapping.parquet",
        "reference.csv": "data/raw/cyp-challenge-TEST-BLINDED.csv",
    }
    spec = cfg["variants"][cfg["variant"]]
    purpose = (
        f"Package D6 {spec['role']} candidate {cfg['variant']} ({spec['formula']}, top "
        f"{cfg['d6_positive_calls']}) from the deployed Sep 15 blind score and run "
        f"{cfg['score_run']} blind scores, with the exact best-scored A4 249-call vector; "
        "no fitting, no upload"
    )
    with tempfile.NamedTemporaryFile("w", suffix=f"-{cfg['variant']}.json", delete=False) as stream:
        json.dump(cfg, stream, indent=2)
        stream.write("\n")
        resolved = Path(stream.name)
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        purpose,
        "--config",
        str(resolved),
        "--input-mode",
        "hardlink",
    ]
    for alias, relative in inputs.items():
        path = root / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        command.extend(["--input", f"{alias}={path}"])
    try:
        subprocess.run(
            command + ["--", "python", "scripts/package_cyp_d6_gate_candidate.py"],
            cwd=root,
            check=True,
        )
    finally:
        resolved.unlink()


if __name__ == "__main__":
    main()
