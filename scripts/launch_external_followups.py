"""Sequential GPU follow-ups after the initial four-objective comparison."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
shared = root / "artifacts/runs/shared"
nested = shared / "20260912T052549384296Z-bac5fb3a/outputs/nested-log2fc.parquet"
classifier = shared / "20260912T051543210441Z-a3ea2ec3/outputs"
trees = shared / "20260912T053008533415Z-45e5a1eb/outputs/external-scores.parquet"
commands = [
    [
        "scripts/launch_external_readout.py",
        "--features",
        str(classifier / "external-embeddings.parquet"),
        "--method",
        "external3_class_emb_logfc_tabpfn",
        "--backend",
        "tabpfn",
        "--nested",
        str(nested),
    ],
    [
        "scripts/launch_external_readout.py",
        "--features",
        str(trees),
        "--method",
        "external_tree_logfc_tabpfn",
        "--backend",
        "tabpfn",
        "--prefix",
        "score_",
        "--nested",
        str(nested),
    ],
    [
        "scripts/launch_external_finetune.py",
        "--checkpoint",
        str(classifier / "best.pt"),
        "--tag",
        "external3_class_mt_v1",
    ],
    [
        "scripts/launch_external_program.py",
        "--dataset-run",
        "20260912T051513413487Z-ac9df6b5",
        "--tag",
        "external3_cond_v1",
        "--modes",
        "classification",
        "mixed",
        "--conditioned",
        "--scaling",
        "pooled_shrinkage",
    ],
]
for command in commands:
    print("START", " ".join(command), flush=True)
    subprocess.run(["pixi", "run", "python", *command], cwd=root, check=True)
