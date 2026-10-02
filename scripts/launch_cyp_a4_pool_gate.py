"""Launch the A4 pool and activity-gate variant evaluation against the equal-seven ensemble."""

import argparse
import subprocess
from pathlib import Path

DATA = "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
SPLITS = "artifacts/runs/shared/20260909T135130268493Z-edba6fc0/outputs"
CHEMELEON = "artifacts/runs/shared/20260910T012630927160Z-e9576045/outputs"
EQUAL7_OOF = "artifacts/runs/classification/20260926T150746748801Z-e74748cd/outputs"
EQUAL7_BLIND = "artifacts/runs/classification/20260915T125939829893Z-1d921ffc/outputs"
SUBS = "artifacts/submissions/classification"
INPUTS = {
    "compounds.parquet": f"{DATA}/compounds.parquet",
    "measurements.parquet": f"{DATA}/measurements.parquet",
    "aliases.parquet": f"{DATA}/aliases.parquet",
    "test_mapping.parquet": f"{DATA}/test_mapping.parquet",
    "assignments.parquet": f"{SPLITS}/assignments.parquet",
    "primary-screen.csv": "data/raw/cyp-challenge-single-concentration-TRAIN.csv",
    "chemeleon-embeddings.parquet": f"{CHEMELEON}/embeddings.parquet",
    "equal7-all-pool-oof.parquet": f"{EQUAL7_OOF}/all-pool-oof.parquet",
    "equal7-blind-scores.parquet": f"{EQUAL7_BLIND}/test-scores.parquet",
    "a4-payload-sep12.parquet": f"{SUBS}/20260912T122715261670Z-363bf7b8/predictions.parquet",
    "a4-payload-sep13.parquet": f"{SUBS}/20260913T113807374506Z-77f4825a/predictions.parquet",
    "a4-payload-sep14.parquet": f"{SUBS}/20260914T020759952207Z-37dd7bd8/predictions.parquet",
    "a4-payload-sep15.parquet": f"{SUBS}/20260915T130004513006Z-75b25070/predictions.parquet",
    "sep15-submission.parquet": f"{SUBS}/20261002T092348406577Z-4ea606c1/predictions.parquet",
}


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/classification/a4-pool-gate-v1.json")
    args = parser.parse_args()
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        "A4 TDI: fit official-label, pseudo-negative pool, activity-gated and CheMeleon "
        "models on the fixed umap/scaffold folds, evaluate equal-seven combinations and "
        "gate filters against the submitted ensemble on official labels with nested rate "
        "selection and paired group bootstrap, and check the two historical A4 "
        "transitions for concordance; no upload",
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
        command + ["--", "python", "scripts/evaluate_cyp_a4_pool_gate.py"],
        cwd=root,
        check=True,
    )


if __name__ == "__main__":
    main()
