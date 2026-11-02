"""Launch the frozen-embedding add test for the D6 TDI score, one embedding family per run."""

import subprocess
import sys
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
EMBEDDINGS = {
    "monroe": "shared/20261029T045000189464Z-11b41494/outputs/external-embeddings.parquet",
    "minimol": "shared/20261029T044801948790Z-dd4df92f/outputs/external-embeddings.parquet",
    "orbmol": "regression/20261029T124544572615Z-c6b401f4/inputs/embeddings.parquet",
    "monroe-pretrained": "regression/20261029T115010020552Z-7bae1535/inputs/embeddings.parquet",
}


def launch(token):
    label, _, variant = token.partition(":")
    config = f"d6-embedding-add-{label}{'-' + variant if variant else ''}-v1.json"
    inputs = {
        "testlike-pool-oof.parquet": runs
        / "classification/20261002T215345688167Z-7920b49a/outputs/oof_scores.parquet",
        "assignments.parquet": runs
        / "classification/20261002T221918697661Z-2d8d20f5/inputs/assignments.parquet",
        "analog-pairs.parquet": root
        / "artifacts/tracked/d6-analog-expansion-validation-20260924/anchor-analog-pairs.parquet",
        "embeddings.parquet": runs / EMBEDDINGS[label],
    }
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        f"Add frozen {label} vectors to the deployed L4 D6 score; paired bootstrap comparison",
        "--config",
        str(root / "configs/classification" / config),
    ]
    for alias, path in inputs.items():
        command += ["--input", f"{alias}={path}"]
    subprocess.run(
        command + ["--", "python", "scripts/evaluate_cyp_d6_cypmol_add.py"], cwd=root, check=True
    )


if __name__ == "__main__":
    for name in sys.argv[1:] or list(EMBEDDINGS):
        launch(name)
