"""Launch the A4 equal-seven OOF stratification and blind-population reweighting."""

import argparse
import subprocess
from pathlib import Path

DATA = "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
SPLITS = "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs"
EQUAL7_OOF = "artifacts/runs/classification/20260926T150746748801Z-e74748cd/outputs"
EQUAL7_BLIND = "artifacts/runs/classification/20260915T125939829893Z-1d921ffc/outputs"
SUBMISSION = "artifacts/submissions/classification/20261002T092348406577Z-4ea606c1"
INPUTS = {
    "compounds.parquet": f"{DATA}/compounds.parquet",
    "measurements.parquet": f"{DATA}/measurements.parquet",
    "aliases.parquet": f"{DATA}/aliases.parquet",
    "test_mapping.parquet": f"{DATA}/test_mapping.parquet",
    "assignments.parquet": f"{SPLITS}/assignments.parquet",
    "primary-screen.csv": "data/raw/cyp-challenge-single-concentration-TRAIN.csv",
    "equal7-all-pool-oof.parquet": f"{EQUAL7_OOF}/all-pool-oof.parquet",
    "equal7-blind-scores.parquet": f"{EQUAL7_BLIND}/test-scores.parquet",
    "sep15-submission.parquet": f"{SUBMISSION}/predictions.parquet",
}


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/classification/a4-oof-strata-v1.json")
    args = parser.parse_args()
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        "A4 TDI: stratify the submitted equal-seven complete-pool UMAP OOF by primary "
        "log2fc bin, label provenance and nearest-blind Tanimoto quartile, then reweight "
        "to the blind population with nearest-neighbour transfer and a fold-excluded "
        "log2fc gate; report reweighted prevalence, AUROC, AP and the 200-350 call curve",
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
        command + ["--", "python", "scripts/evaluate_cyp_a4_oof_strata.py"],
        cwd=root,
        check=True,
    )


if __name__ == "__main__":
    main()
