"""Bind completed matching runs to an immutable audit configuration and execute it."""

import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main():
    base = json.loads((ROOT / "configs/external/training-v1.json").read_text())
    selected = {}
    for folder in sorted((ROOT / "artifacts/runs/regression").iterdir()):
        if not (folder / "outputs/summary.json").exists() or not (folder / "result.json").exists():
            continue
        summary = json.loads((folder / "outputs/summary.json").read_text())
        cfg = summary.get("config", {})
        if any(cfg.get(k) != v for k, v in base.items()):
            continue
        if json.loads((folder / "result.json").read_text())["status"] != "completed":
            continue
        key = cfg["variant"], cfg["fold"]
        if key in selected:
            raise ValueError(f"Ambiguous completed run for {key}")
        selected[key] = folder
    assert set(selected) == {(v, f) for v in ["control", "external"] for f in range(5)}
    files = []

    def add(path):
        name = f"{len(files):03d}-{path.name}"
        files.append(path)
        return name

    shared = ROOT / "artifacts/runs/shared"
    cfg = dict(
        labels=add(shared / base["dataset_run"] / "outputs/measurements.parquet"),
        splits=add(shared / base["split_run"] / "outputs/assignments.parquet"),
        scorer=add(ROOT / "artifacts/metric-reference/custom_scoring_functions.py"),
        test_reference=add(
            shared / "20260910T132141471698Z-343eb5ec/outputs/ensemble_probabilities.parquet"
        ),
        bootstrap_seed=20260911,
        bootstrap_replicates=2000,
        variants={
            "baseline": [
                {"oof": add(shared / "20260910T123856868493Z-243a1e23/outputs/oof.parquet")}
            ]
        },
    )
    for variant in ["control", "external"]:
        cfg["variants"][variant] = []
        for fold in range(5):
            folder = selected[(variant, fold)]
            cfg["variants"][variant].append(
                dict(
                    run_id=folder.name,
                    fold=fold,
                    oof=add(folder / "outputs/oof.parquet"),
                    summary=add(folder / "outputs/summary.json"),
                    membership=add(folder / "outputs/membership.parquet"),
                    test=add(folder / "outputs/test-predictions.parquet"),
                )
            )
    path = ROOT / "configs/external/comparison-audit-v1.json"
    with path.open("x") as stream:
        json.dump(cfg, stream, indent=2)
        stream.write("\n")
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        "Audit matched five-fold external-data comparison and paired UMAP-cluster uncertainty",
        "--config",
        str(path),
    ]
    for file in files:
        command += ["--input", str(file)]
    subprocess.run(
        command + ["--", "python", "scripts/audit_external_comparison.py"], cwd=ROOT, check=True
    )


if __name__ == "__main__":
    main()
