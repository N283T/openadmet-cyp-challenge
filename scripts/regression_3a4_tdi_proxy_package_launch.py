"""Launch the fixed-weight 3A4 TDI-proxy member swap evaluation and payload build."""

import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
bank = runs / "regression/20260917T130718808881Z-dfa02730/inputs"
inputs = {path.name: path for path in bank.iterdir() if path.suffix == ".parquet"}
inputs.pop("current-control.parquet", None)
inputs.update(
    {
        "custom_scoring_functions.py": bank / "custom_scoring_functions.py",
        "control-3a4.parquet": runs
        / "regression/20260917T130718808881Z-dfa02730/outputs/oof.parquet",
        "parent-assignments.parquet": runs
        / "regression/20260911T104703655337Z-a3923255/outputs/parent-assignments.parquet",
        "analog-pairs.parquet": runs
        / "regression/20260927T054145210148Z-ba3de34b/inputs/analog-pairs.parquet",
        "proxy-fast-oof.parquet": runs
        / "regression/20261029T053332281624Z-78ab3369/outputs/oof.parquet",
        "proxy-full-oof.parquet": runs
        / "regression/20261029T054151324166Z-7f7fe3da/outputs/oof.parquet",
        "proxy-fast-test.parquet": runs
        / "regression/20261029T054252727488Z-28487958/outputs/test-predictions.parquet",
        "proxy-full-test.parquet": runs
        / "regression/20261029T054340496124Z-4a4c01b0/outputs/test-predictions.parquet",
        "atomic-test.parquet": runs
        / "regression/20260917T131642179677Z-cf8c2eb0/outputs/atomic-test-predictions.parquet",
        "sent-predictions.csv": root
        / "artifacts/submissions/regression/20261001T123022061267Z-ad1689ae/predictions.csv",
        "test-mapping.parquet": runs
        / "shared/20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
    }
)
command = [
    "pixi",
    "run",
    "cyp",
    "run",
    "--track",
    "regression",
    "--purpose",
    (
        "Swap the two anchor59 TabPFN 3A4 members for TDI-proxy refits at fixed deployed "
        "weights; evaluate on full and analog OOF and build the 3A4-only payload change"
    ),
    "--config",
    str(root / "configs/regression/3a4-tdi-proxy-package-v1.json"),
]
for alias, path in sorted(inputs.items()):
    command += ["--input", f"{alias}={path}"]
subprocess.run(
    command + ["--", "python", "scripts/regression_3a4_tdi_proxy_package.py"], cwd=root, check=True
)
