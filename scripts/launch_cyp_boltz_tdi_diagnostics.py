"""Record compact labelled-only TDI integration diagnostics without new fits."""

import json
import subprocess
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    prior = root / "artifacts/runs/classification/20261001T065532999859Z-47297bc3"
    manifest = json.loads((prior / "manifest.json").read_text())
    inputs = {Path(item["snapshot"]).name: Path(item["snapshot"]) for item in manifest["inputs"]}
    inputs.update(
        {
            "placement-oof.parquet": prior / "outputs/oof.parquet",
            "measurements.parquet": root
            / "artifacts/runs/classification/20261001T023655009132Z-91096754"
            / "inputs/measurements.parquet",
            "assignments.parquet": root
            / "artifacts/runs/classification/20261001T023655009132Z-91096754"
            / "inputs/assignments.parquet",
            "tensor-manifest.parquet": root / "artifacts/boltz-trunk/post-r1/manifest.parquet",
        }
    )
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        "Diagnose saved TDI trunk addition, replacement, diversity, "
        "and operating-boundary signals on the labelled pool",
        "--config",
        str(root / "configs/boltz-trunk/tdi-partial-pool-diagnostics-v1.json"),
    ]
    for alias, path in inputs.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        command.extend(["--input", f"{alias}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/evaluate_cyp_boltz_tdi_diagnostics.py"],
        cwd=root,
        check=True,
    )


if __name__ == "__main__":
    main()
