"""Record exact source artifacts for CPU-only historical feature extraction."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
b = root / "artifacts/runs"
old = b / "shared/20260910T132141471698Z-343eb5ec"
new = b / "classification/20260912T132823500863Z-130ac6af"
o = json.loads((old / "config.json").read_text())
n = json.loads((new / "config.json").read_text())
inputs = []
cache = {}


def add(p):
    if p not in cache:
        cache[p] = f"{len(inputs):03d}-{p.name}"
        inputs.append(p)
    return cache[p]


cfg = {k: add(new / "inputs" / n[k]) for k in ["labels", "assignments"]}
cfg["test"] = add(b / "shared/20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet")
for k in ["embedding", "mordred"]:
    cfg[k] = [add(new / "inputs" / p) for p in n[k]]
cfg["auxiliary"] = []
seen = set()
for fold in o["folds"]:
    for s in fold["auxiliary"]:
        meta = json.loads((old / "inputs" / s["summary"]).read_text())
        excluded = tuple(meta["config"]["excluded_folds"])
        if excluded in seen:
            continue
        seen.add(excluded)
        cfg["auxiliary"].append(
            dict(
                excluded=excluded,
                summary=add(old / "inputs" / s["summary"]),
                predictions=[add(old / "inputs" / p) for p in s["predictions"]],
            )
        )
for s in n["auxiliary"]:
    cfg["auxiliary"].append(
        dict(
            excluded=s["excluded_folds"],
            summary=add(new / "inputs" / s["summary"]),
            predictions=[add(new / "inputs" / p) for p in s["predictions"]],
        )
    )
cfg["readouts"] = []
for fold in o["folds"]:
    outer = fold["fold"]
    for ep in ["CYP2D6_is_TDI", "CYP3A4_is_TDI"]:
        paths = fold["endpoints"][ep]
        cfg["readouts"].append(
            dict(
                outer=outer,
                endpoint=ep,
                phase="outer",
                **{
                    k: add(old / "inputs" / paths[v])
                    for k, v in [
                        ("features", "validation_features"),
                        ("pca", "pca"),
                        ("imputer", "imputer"),
                        ("columns", "columns"),
                    ]
                },
            )
        )
        cfg["readouts"].append(
            dict(
                outer=outer,
                endpoint=ep,
                phase="inner",
                preprocessors=add(new / "outputs" / f"{ep}-outer{outer}-preprocessors.joblib"),
            )
        )
p = root / "configs/classification-stack-icl"
p.mkdir(exist_ok=True)
conf = p / "prepare.json"
conf.write_text(json.dumps(cfg, indent=2))
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    (
        "Extract exact historical stack outer inner and test matrices for "
        "TabICL classifier-only comparison"
    ),
    "--config",
    str(conf),
]
for p in inputs:
    cmd += ["--input", str(p)]
subprocess.run(
    cmd + ["--", "python", "scripts/classification_stack_icl_prepare.py"], cwd=root, check=True
)
