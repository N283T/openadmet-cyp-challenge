"""Sequential recorded source expansion, auxiliary and source-ablation experiments."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
shared = root / "artifacts/runs/shared"
checkpoint = shared / "20260912T051543210441Z-a3ea2ec3/outputs/best.pt"
emax = root / "artifacts/runs/regression/20260912T055345388545Z-8ab226fa/outputs/targets.parquet"
commands = [
    [
        "scripts/launch_external_readout.py",
        "--features",
        str(shared / "20260912T055819206227Z-c060d197/outputs/external-embeddings.parquet"),
        "--method",
        "jeremy_upstream_medium_tabpfn_diagnostic_v1",
        "--backend",
        "tabpfn",
    ],
    [
        "scripts/launch_external_program.py",
        "--dataset-run",
        "20260912T055348682784Z-be731fe3",
        "--tag",
        "external_az_mixed_v1",
        "--modes",
        "mixed",
        "--scaling",
        "pooled_shrinkage",
        "--finetune",
    ],
    [
        "scripts/launch_external_finetune.py",
        "--checkpoint",
        str(checkpoint),
        "--tag",
        "external3_class_emax_v1",
        "--epochs",
        "160",
        "--patience",
        "20",
        "--extra-targets",
        str(emax),
        "--extra-weight",
        "0.5",
    ],
    [
        "scripts/launch_external_program.py",
        "--dataset-run",
        "20260912T054023889952Z-3ba75913",
        "--tag",
        "pubchem_tox_class_v1",
        "--modes",
        "classification",
    ],
    [
        "scripts/launch_external_program.py",
        "--dataset-run",
        "20260912T055809060996Z-06088efb",
        "--tag",
        "external_adme_medium_v1",
        "--modes",
        "mixed",
        "--architecture",
        "medium",
        "--scaling",
        "pooled_shrinkage",
        "--finetune",
    ],
]
for command in commands:
    print("START", " ".join(command), flush=True)
    subprocess.run(["pixi", "run", "python", *command], cwd=root, check=True)
