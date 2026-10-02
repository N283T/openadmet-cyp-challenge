"""Record failure-linked decision diagnostics and a concrete prospective D6 payload."""

import subprocess
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    runs = root / "artifacts/runs/classification"
    new = runs / "20261002T021327768010Z-c37abd9c/outputs"
    data = root / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
    inputs = {
        f"new-{name}": new / name
        for name in ["oof.parquet", "summary.json", "operating-points.csv", "blind-folds.parquet"]
    }
    inputs.update(
        {
            "control-umap.parquet": runs
            / "20260924T040301085704Z-914bb06c/outputs/all-pool-oof.parquet",
            "control-scaffold.parquet": runs
            / "20260924T041452692950Z-e17778e3/outputs/all-pool-oof.parquet",
            "ap-umap.parquet": runs
            / "20260924T035023873825Z-b719ba35/outputs/all-pool-oof.parquet",
            "ap-scaffold.parquet": runs
            / "20260924T041348966036Z-8db19b4b/outputs/all-pool-oof.parquet",
            "trunk-oof.parquet": root
            / "artifacts/tracked/cyp-boltz-tdi-diagnostics-20261001T133225940252Z-502d403a"
            / "selected-oof.parquet",
            "live-summary.json": runs / "20261001T143048820695Z-05fb009d/outputs/summary.json",
            "last-sent.parquet": root
            / "artifacts/submissions/classification/20261001T141041310215Z-c7faa79f"
            / "predictions.parquet",
            "best-control.parquet": root
            / "artifacts/submissions/classification/20260924T050705396203Z-10e1d3ad"
            / "predictions.parquet",
            "reference.csv": root / "data/raw/cyp-challenge-TEST-BLINDED.csv",
            "aliases.parquet": data / "aliases.parquet",
            "historical-replay.md": root / "docs/analysis/2026-09-24-d6-ap-full-pool-replay.md",
        }
    )
    for path in new.glob("blind-*.parquet"):
        if path.name != "blind-folds.parquet":
            inputs[path.name] = path
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        "Link D6 AP/Live failures to model, full-pool decision and inference effects; "
        "freeze independent probability-policy submission probe for review",
        "--config",
        str(root / "configs/classification/d6-causal-probe-review-v1.json"),
    ]
    for alias, path in inputs.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        command.extend(["--input", f"{alias}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/review_cyp_d6_causal_probe.py"], cwd=root, check=True
    )


if __name__ == "__main__":
    main()
