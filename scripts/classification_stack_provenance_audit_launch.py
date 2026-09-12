"""Snapshot final reconstruction evidence for independent non-training audit."""

import argparse
import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
a = argparse.ArgumentParser()
a.add_argument("--run", required=True)
args = a.parse_args()
run = runs / "classification" / args.run
src = json.loads((run / "config.json").read_text())
inputs = []


def add(p):
    assert p.is_file(), p
    inputs.append(p)
    return len(inputs) - 1


cfg = {
    k: add(run / "inputs" / src[k])
    for k in ["labels", "assignments", "old_outer", "old_test", "old_mean"]
}
for k in ["embedding", "mordred"]:
    cfg[k] = [add(run / "inputs" / p) for p in src[k]]
cfg["targets"] = add(
    next((runs / "shared" / src["auxiliary"][0]["run"] / "inputs").glob("*targets.parquet"))
)
cfg["auxiliary"] = []
for s in src["auxiliary"]:
    d = runs / "shared" / s["run"] / "outputs"
    item = dict(run=s["run"], excluded_folds=s["excluded_folds"])
    for key, file in [
        ("summary", "summary.json"),
        ("metadata", "metadata.json"),
        ("history", "history.jsonl"),
        ("membership", "membership.parquet"),
        ("checkpoint", "best.pt"),
    ]:
        item[key] = add(d / file)
    item["predictions"] = [add(run / "inputs" / p) for p in s["predictions"]]
    cfg["auxiliary"].append(item)
for k, file in [
    ("inner", "inner-probabilities.parquet"),
    ("metrics", "inner-metrics.csv"),
    ("summary", "summary.json"),
    ("new_outer", "outer-oof.parquet"),
    ("new_test", "fold-test-probabilities.parquet"),
    ("new_mean", "test-probabilities.parquet"),
]:
    cfg[k] = add(run / "outputs" / file)
cfg["readouts"] = []
for outer in range(5):
    for ep in ["CYP2D6_is_TDI", "CYP3A4_is_TDI"]:
        stem = f"{ep}-outer{outer}"
        cfg["readouts"].append(
            dict(
                outer_fold=outer,
                endpoint=ep,
                membership=add(run / "outputs" / f"{stem}-membership.parquet"),
                preprocessor=add(run / "outputs" / f"{stem}-preprocessors.joblib"),
            )
        )
config = root / "configs/classification-stack-inner/provenance-audit.json"
config.write_text(json.dumps(cfg, indent=2))
cmd = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "classification",
    "--purpose",
    "Independent historical TabPFN nested-inner auxiliary/PCA/descriptor/identity provenance audit",
    "--config",
    str(config),
]
for p in inputs:
    cmd += ["--input", str(p)]
subprocess.run(
    cmd + ["--", "python", "scripts/classification_stack_provenance_audit.py"], cwd=root, check=True
)
