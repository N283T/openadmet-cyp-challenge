"""Refit the two deployed anchor59 TabPFN 3.5 3A4 members with TDI-proxy labels.

Every input except measurements is the exact snapshot of the original September 17
full-data refit; the configuration only renames the method and restricts it to 3A4.
"""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
measurements = runs / "shared/20261029T054205583446Z-3f976428/outputs/measurements.parquet"
originals = {
    "fast": "regression/20260917T124004117613Z-279d6c0f",
    "full": "regression/20260917T124342046803Z-4f691c7d",
}
for variant, source in originals.items():
    base = runs / source
    cfg = json.loads((base / "config.json").read_text())
    cfg["method"] = cfg["method"].replace("_fullfit_v1", "_tdi_proxy_fullfit_v1")
    cfg["endpoints"] = ["CYP3A4_pIC50_direct_inhibition"]
    cfg["expected_label_rows"] = 6525 + 1249
    config = root / f"configs/regression/3a4-tdi-proxy-fullfit-{variant}-v1.json"
    config.write_text(json.dumps(cfg, indent=2) + "\n")
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "regression",
        "--purpose",
        f"Full-data 3A4 refit of {cfg['method']} with TDI-proxy direct labels",
        "--config",
        str(config),
    ]
    for path in sorted((base / "inputs").iterdir()):
        if path.is_file():
            source_path = measurements if path.name == "measurements.parquet" else path
            command += ["--input", f"{path.name}={source_path}"]
    subprocess.run(
        command + ["--", "python", "scripts/full_submission_readout.py"], cwd=root, check=True
    )
