"""Record full-OOF raw LightGBM audits for downstream-guided checkpoints."""

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OBJECTIVE = ROOT / "artifacts/runs/regression/20260916T120842454770Z-b5a003ef/outputs"
SOURCE = ROOT / "artifacts/runs/regression/20260916T123149468017Z-68f542a6/outputs"
CONFIG = ROOT / "configs/owned-encoder-v2/downstream_task_design_full_oof_v1.json"
MEASUREMENTS = (
    ROOT / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs" / "measurements.parquet"
)
ASSIGNMENTS = (
    ROOT
    / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
    / "parent-assignments.parquet"
)

files = {
    "objective-source-uniform.parquet": OBJECTIVE
    / "trial-00-control_source_uniform_mixed/external-embeddings.parquet",
    "objective-task-uniform.parquet": OBJECTIVE
    / "trial-01-task_uniform_mixed/external-embeddings.parquet",
    "objective-sqrt-mixed.parquet": OBJECTIVE / "trial-02-sqrt_mixed/external-embeddings.parquet",
    "objective-proportional-mixed.parquet": OBJECTIVE
    / "trial-03-proportional_mixed/external-embeddings.parquet",
    "objective-sqrt-weak-aux.parquet": OBJECTIVE
    / "trial-04-sqrt_weak_aux/external-embeddings.parquet",
    "objective-sqrt-weak-classification.parquet": OBJECTIVE
    / "trial-05-sqrt_weak_classification/external-embeddings.parquet",
    "objective-sqrt-no-classification.parquet": OBJECTIVE
    / "trial-06-sqrt_no_classification/external-embeddings.parquet",
    "source-no-tdc.parquet": SOURCE / "trial-00-sqrt_no_tdc/external-embeddings.parquet",
    "source-no-pubchem.parquet": SOURCE / "trial-01-sqrt_no_pubchem/external-embeddings.parquet",
    "source-no-tox21.parquet": SOURCE / "trial-02-sqrt_no_tox21/external-embeddings.parquet",
    "source-no-polaris.parquet": SOURCE / "trial-03-sqrt_no_polaris/external-embeddings.parquet",
    "source-no-expansionrx.parquet": SOURCE
    / "trial-04-sqrt_no_expansionrx/external-embeddings.parquet",
    "source-no-pxr.parquet": SOURCE / "trial-05-sqrt_no_pxr/external-embeddings.parquet",
    "source-cyp-only.parquet": SOURCE
    / "trial-06-sqrt_cyp_sources_only/external-embeddings.parquet",
}

command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Full-OOF target-specific audit of downstream-guided task designs",
    "--config",
    str(CONFIG),
    "--input",
    f"measurements.parquet={MEASUREMENTS}",
    "--input",
    f"parent-assignments.parquet={ASSIGNMENTS}",
    "--input",
    "custom_scoring_functions.py="
    + str(ROOT / "artifacts/metric-reference/custom_scoring_functions.py"),
]
for alias, path in files.items():
    command.extend(["--input", f"{alias}={path}"])
command.extend(
    [
        "--input",
        "previous-owned.parquet="
        + str(
            ROOT
            / "artifacts/runs/regression/20260916T085659227651Z-a4bc5a2c/outputs"
            / "external-oof.parquet"
        ),
        "--input",
        "jeremy.parquet="
        + str(
            ROOT
            / "artifacts/runs/regression/20260912T055849013223Z-f16a5dfd/outputs"
            / "external-oof.parquet"
        ),
        "--",
        "python",
        "scripts/focused_cyp14_architecture_lgbm_screen.py",
    ]
)
subprocess.run(command, cwd=ROOT, check=True)
