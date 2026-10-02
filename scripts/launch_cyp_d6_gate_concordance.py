"""Launch the D6 activity-gate historical concordance and legacy x gate OOF evaluation."""

import argparse
import subprocess
from pathlib import Path

DATA = "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
SPLITS = "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs"
POOL_RUN = "artifacts/runs/classification/20261002T215345688167Z-7920b49a/outputs"
SEP15_BLIND = "artifacts/runs/classification/20260924T045344348206Z-7e2db311/inputs"
SUBMISSIONS = "artifacts/submissions/classification"
PAYLOADS = {
    "payload_sep13.parquet": "20260913T113807374506Z-77f4825a",
    "payload_sep14.parquet": "20260914T020759952207Z-37dd7bd8",
    "payload_sep15.parquet": "20260915T130004513006Z-75b25070",
    "payload_sep18.parquet": "20260918T094451131332Z-f23dd61b",
    "payload_sep23.parquet": "20260923T180233070503Z-552d9af4",
    "payload_oct01.parquet": "20261001T141041310215Z-c7faa79f",
    "payload_oct02.parquet": "20261002T092348406577Z-4ea606c1",
}
INPUTS = {
    "oof_scores.parquet": f"{POOL_RUN}/oof_scores.parquet",
    "test_scores.parquet": f"{POOL_RUN}/test_scores.parquet",
    "sep15_blind_scores.parquet": f"{SEP15_BLIND}/blind-scores.parquet",
    "aliases.parquet": f"{DATA}/aliases.parquet",
    "test_mapping.parquet": f"{DATA}/test_mapping.parquet",
    "assignments.parquet": f"{SPLITS}/assignments.parquet",
    "primary-screen.csv": "data/raw/cyp-challenge-single-concentration-TRAIN.csv",
    **{alias: f"{SUBMISSIONS}/{run}/predictions.parquet" for alias, run in PAYLOADS.items()},
}
PURPOSE = (
    "D6 activity gate historical LB concordance: permutation test of removed versus added "
    "calls across the seven submitted D6 payloads, legacy x gate OOF combinations against the "
    "legacy Sep 15 refit control on umap/scaffold, blind swap inventory and a strict "
    "pseudo-negative sensitivity; read-only, no candidate or upload"
)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/classification/d6-gate-concordance-v1.json")
    args = parser.parse_args()
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
    for alias, relative in INPUTS.items():
        path = root / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        command.extend(["--input", f"{alias}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/evaluate_cyp_d6_gate_concordance.py"],
        cwd=root,
        check=True,
    )


if __name__ == "__main__":
    main()
