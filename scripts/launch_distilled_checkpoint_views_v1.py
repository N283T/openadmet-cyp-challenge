"""Export multiple auditable feature views from all distillation checkpoints."""

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "artifacts/runs/regression"
CHALLENGE = ROOT / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
DISTILLED = RUNS / "20260917T051118981633Z-a9758b6c/outputs"
PUBCHEM = RUNS / "20260917T091542338422Z-b0467aeb/outputs"
inputs = {
    "compounds.parquet": CHALLENGE / "compounds.parquet",
    "scratch-teacher-only.pt": DISTILLED / "trial-00-scratch_teacher_only/checkpoint.pt",
    "owned-init-hybrid.pt": DISTILLED / "trial-01-owned_init_hybrid/checkpoint.pt",
    "jeremy-init-task-anchor.pt": DISTILLED / "trial-02-jeremy_init_task_anchor/checkpoint.pt",
    "jeremy-pubchem-only.pt": PUBCHEM / "trial-00-jeremy_pubchem_only/checkpoint.pt",
    "jeremy-pubchem-anchor010.pt": PUBCHEM
    / "trial-01-jeremy_pubchem_anchor010/checkpoint.pt",
    "jeremy-pubchem-anchor030.pt": PUBCHEM
    / "trial-02-jeremy_pubchem_anchor030/checkpoint.pt",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    "Export raw, final-hidden, and external-task predictor views from owned checkpoints",
    "--config",
    str(ROOT / "configs/external-training/distilled-checkpoint-views-v1.json"),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/export_distilled_checkpoint_views.py"],
    cwd=ROOT,
    check=True,
)
