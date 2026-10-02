"""Record label-selection audit against successful and failed historical changes."""

import subprocess
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    runs = root / "artifacts/runs"
    shared = runs / "shared/20260909T132734119233Z-a7063ccd/outputs"
    inputs = {
        name: shared / name
        for name in ["compounds.parquet", "aliases.parquet", "measurements.parquet"]
    }
    inputs["assignments.parquet"] = (
        runs / "shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet"
    )
    inputs["primary-screen.csv"] = root / "data/raw/cyp-challenge-single-concentration-TRAIN.csv"
    for alias, run, filename in [
        ("diverse-oof.parquet", "20260913T112957931554Z-32d00f6e", "oof.parquet"),
        ("fixed-oof.parquet", "20260913T123045515718Z-0df22d54", "oof.parquet"),
        ("control-umap.parquet", "20260924T040301085704Z-914bb06c", "all-pool-oof.parquet"),
        ("control-scaffold.parquet", "20260924T041452692950Z-e17778e3", "all-pool-oof.parquet"),
        ("ap-umap.parquet", "20260924T035023873825Z-b719ba35", "all-pool-oof.parquet"),
        ("ap-scaffold.parquet", "20260924T041348966036Z-8db19b4b", "all-pool-oof.parquet"),
    ]:
        inputs[alias] = runs / "classification" / run / "outputs" / filename
    inputs["official-assay.md"] = (
        root / "references/official/cyp/assay-dataset/20260907T135802Z/snapshot.md"
    )
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        "Test primary-conditioned class-label selection against successful and failed "
        "historical D6 changes without assigning unknown labels",
        "--config",
        str(root / "configs/classification/d6-primary-selection-audit-v1.json"),
    ]
    for alias, path in inputs.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        command.extend(["--input", f"{alias}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/audit_cyp_d6_primary_selection.py"],
        cwd=root,
        check=True,
    )


if __name__ == "__main__":
    main()
