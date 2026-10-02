"""Record full-development recipe selection for an unsubmitted D6 call review."""

import json
import subprocess
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    source = root / "artifacts/runs/classification/20261001T133225940252Z-502d403a"
    cfg = json.loads((source / "config.json").read_text())
    cfg.update({"diagnostic_source": source.name, "positive_calls": 49})
    config = root / "configs/boltz-trunk/d6-ap-probe-selection-v1.json"
    if config.exists():
        if json.loads(config.read_text()) != cfg:
            raise ValueError("Existing configuration differs; create a new version")
    else:
        with config.open("x") as stream:
            json.dump(cfg, stream, indent=2)
            stream.write("\n")
    manifest = json.loads((source / "manifest.json").read_text())
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        "Select diagnosed D6 AP-block recipe for review of 750 blind calls; no upload",
        "--config",
        str(config),
    ]
    for item in manifest["inputs"]:
        path = Path(item["snapshot"])
        command.extend(["--input", f"{path.name}={path}"])
    subprocess.run(
        command + ["--", "python", "scripts/select_cyp_boltz_d6_probe.py"], cwd=root, check=True
    )


if __name__ == "__main__":
    main()
