"""Launch all-training fits and package three D6 checkpoint candidates."""

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "artifacts/runs/regression"
SHARED = ROOT / "artifacts/runs/shared"
VIEWS = RUNS / "20260917T111518065286Z-801ab060/outputs"
BROAD = RUNS / "20260917T111705603506Z-e23f10b5/outputs"
DIVERSE = RUNS / "20260917T112238487125Z-7178add0/outputs"
inputs = {
    "measurements.parquet": SHARED
    / "20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    "aliases.parquet": SHARED / "20260909T132734119233Z-a7063ccd/outputs/aliases.parquet",
    "reference.csv": ROOT / "data/raw/cyp-challenge-TEST-BLINDED.csv",
    "base-candidate.csv": ROOT
    / "artifacts/candidates/regression/20260914T024735834387Z-f87e9914/predictions.csv",
    "custom_scoring_functions.py": ROOT
    / "artifacts/metric-reference/custom_scoring_functions.py",
    "jeremy-frozen-test.parquet": RUNS
    / "20260912T074409483212Z-f3024998/outputs/test-predictions.parquet",
    "jeremy-scores-test.parquet": RUNS
    / "20260912T074557646349Z-da857069/outputs/test-predictions.parquet",
    "jeremy-finetuned-test.parquet": RUNS
    / "20260912T074624192288Z-bd58b762/outputs/test-predictions.parquet",
    "chemeleon-raw.parquet": SHARED
    / "20260912T061231617816Z-299439e8/outputs/external-embeddings.parquet",
    "full-aux.parquet": SHARED
    / "20260912T074339967571Z-ff194fcf/outputs/full-aux.parquet",
    "tabicl-regressor.ckpt": SHARED
    / "20260911T093102177307Z-beddd36f/outputs/tabicl-regressor-v2-20260212.ckpt",
    "owned-raw-scores.parquet": VIEWS / "owned_init_hybrid-raw-scores.parquet",
    "anchor59-raw-scores.parquet": VIEWS / "jeremy_init_task_anchor-raw-scores.parquet",
    "anchor59-scores.parquet": VIEWS / "jeremy_init_task_anchor-scores.parquet",
    "pubonly-raw.parquet": VIEWS / "jeremy_pubchem_only-raw.parquet",
    "scratch-raw.parquet": VIEWS / "scratch_teacher_only-raw.parquet",
    "oof-jeremy-frozen.parquet": RUNS
    / "20260912T060059705105Z-2f0fbfbd/outputs/external-oof.parquet",
    "oof-jeremy-scores.parquet": RUNS
    / "20260912T072303592078Z-07f2749a/outputs/external-oof.parquet",
    "oof-jeremy-finetuned.parquet": RUNS
    / "20260912T072056445264Z-e6b6cbcb/outputs/external-oof.parquet",
    "oof-chemeleon-raw.parquet": RUNS
    / "20260916T060715739725Z-b0ddf3c8/outputs/external-oof.parquet",
    "oof-owned-raw-scores.parquet": BROAD
    / "owned_init_hybrid_raw_scores__lgbm-oof.parquet",
    "oof-anchor59-raw-scores.parquet": BROAD
    / "jeremy_init_task_anchor_raw_scores__lgbm-oof.parquet",
    "oof-anchor59-scores.parquet": BROAD
    / "jeremy_init_task_anchor_scores__lgbm-oof.parquet",
    "oof-pubonly-raw.parquet": BROAD / "jeremy_pubchem_only_raw__lgbm-oof.parquet",
    "oof-scratch-catboost.parquet": DIVERSE
    / "scratch_teacher_only_raw__catboost-oof.parquet",
    "oof-anchor59-catboost.parquet": DIVERSE
    / "jeremy_init_task_anchor_raw_scores__catboost-oof.parquet",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Full-fit D6 checkpoint readouts and package stable5, checkpoint6 and equal10 candidates",
    "--config",
    str(ROOT / "configs/regression/distilled-d6-fullfit-candidates-v1.json"),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/fullfit_distilled_d6_candidates.py"],
    cwd=ROOT,
    check=True,
)
