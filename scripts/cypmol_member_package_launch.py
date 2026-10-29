"""Launch the CYPMol-ensemble member refit and payload build."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
regression, shared = runs / "regression", runs / "shared"
inputs = {
    "embeddings.parquet": shared / "20261029T142347849581Z-2c6a74f4/outputs/external-embeddings.parquet",
    "screen-blind-mean.parquet": regression
    / "20261029T145728905917Z-9d7ea648/outputs/predictions/blind-mean.parquet",
    "base-predictions.csv": regression / "20261029T054536454025Z-220f26c8/outputs/predictions.csv",
    "test-mapping.parquet": shared / "20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
    "tabpfn.ckpt": shared
    / "20260910T122343854890Z-3a03f908/inputs/056-tabpfn-v2.6-regressor-v2.6_default.ckpt",
    "control-1a2.parquet": regression
    / "20260918T034159381060Z-46314b50/outputs/1A2/add_dense_and_fast_pool/oof.parquet",
    "control-2d6.parquet": regression / "20260926T085959697226Z-0dbf57d4/outputs/oof.parquet",
}
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    (
        "Refit CYPMol-ensemble readouts on all labelled 1A2 and D6 rows and blend them at fixed "
        "weights 0.2 and 0.1 into the 3A4 TDI-proxy payload"
    ),
    "--config",
    str(root / "configs/regression/cypmol-member-package-v1.json"),
]
for alias, path in inputs.items():
    command += ["--input", f"{alias}={path}"]
subprocess.run(
    command + ["--", "python", "scripts/cypmol_member_package.py"], cwd=root, check=True
)
