"""Record immutable official snapshots for the scored AP-trunk probe audit."""

import subprocess
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    api = root / "artifacts/api/classification"
    inputs = {
        "receipt.json": root
        / "artifacts/submission-attempts/classification/20261001T141324626858Z-cd9cca73"
        / "result.json",
        "submitted.parquet": root
        / "artifacts/submissions/classification/20261001T141041310215Z-c7faa79f"
        / "predictions.parquet",
        "last.parquet": root
        / "artifacts/submissions/classification/20260926T150948735129Z-6fc9d05b"
        / "predictions.parquet",
        "best.parquet": root
        / "artifacts/submissions/classification/20260924T050705396203Z-10e1d3ad"
        / "predictions.parquet",
    }
    snapshots = {
        "MA-before.csv": "20261001T141311749379Z-8d98e149",
        "D6-before.csv": "20261001T141316010733Z-6e8d9b1a",
        "A4-before.csv": "20261001T141320306325Z-da763930",
        "MA-after.csv": "20261001T142809835192Z-c0266018",
        "D6-after.csv": "20261001T142809834162Z-478c049d",
        "A4-after.csv": "20261001T142809834162Z-9ada1e4e",
        "A4-best.csv": "20260926T151004586043Z-57198edf",
    }
    for alias, identifier in snapshots.items():
        paths = list((api / identifier).rglob("*.csv"))
        if len(paths) != 1:
            raise ValueError(f"Ambiguous snapshot {identifier}")
        inputs[alias] = paths[0]
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        "Audit scored D6 AP-trunk Live loss and exact A4 best-control restoration",
        "--config",
        str(root / "configs/boltz-trunk/classification-live-audit-20261001.json"),
    ]
    for alias, path in inputs.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        command.extend(["--input", f"{alias}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/audit_cyp_boltz_classification_live.py"],
        cwd=root,
        check=True,
    )


if __name__ == "__main__":
    main()
