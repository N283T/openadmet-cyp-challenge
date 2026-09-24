"""Launch full-fit AP mass90 CYP2D6 candidate composition."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
config = root / "configs/classification/d6-ap-mass90-fullfit-v1.json"
fixed = runs / "20260913T123045515718Z-0df22d54"
inputs = {
    "legacy-bank.parquet": runs / "20260915T121009912691Z-73198c15/outputs/prediction-bank.parquet",
    "lgbm-oof.parquet": runs / "20260919T155849130417Z-0394b600/outputs/oof.parquet",
    "lgbm-test.parquet": runs / "20260919T155849130417Z-0394b600/outputs/test.parquet",
    "readout-oof.parquet": runs / "20260923T161449619218Z-a5f9f7ea/outputs/oof.parquet",
    "readout-test.parquet": runs / "20260923T161449619218Z-a5f9f7ea/outputs/test.parquet",
    "tabpfn-oof.parquet": runs / "20260923T161718040624Z-8a9676e7/outputs/oof.parquet",
    "tabpfn-test.parquet": runs / "20260923T161718040624Z-8a9676e7/outputs/test.parquet",
    "mace-test.parquet": runs / "20260914T051014577658Z-7b1b64c2/outputs/test.parquet",
    "activity-test.parquet": runs / "20260915T065610180920Z-f532943e/outputs/test.parquet",
    "mechanism-test.parquet": runs / "20260914T063614433212Z-c0dc3f32/outputs/test.parquet",
    "pfn-test.parquet": runs
    / "20260912T121146867140Z-b0a05986/inputs/005-test-fold-probabilities.parquet",
    "pool-definition.json": root / "configs/classification/d6-flat-common-pool-caruana-v1.json",
    "outer-selections.parquet": runs
    / "20260923T165552531859Z-d49f8c8b/outputs/outer-selections.parquet",
    "adaptive-oof.parquet": runs / "20260923T171350543926Z-9081bf31/outputs/oof.parquet",
    "control-oof.parquet": runs / "20260918T085712373599Z-8829d528/outputs/oof.parquet",
    "control-test.parquet": runs / "20260918T085712373599Z-8829d528/outputs/test.parquet",
    "compounds.parquet": root
    / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    "test-mapping.parquet": fixed / "inputs/002-test_mapping.parquet",
}
for fold in range(5):
    inputs[f"complement-fold{fold}.joblib"] = fixed / f"outputs/complement_nb_a1-fold{fold}.joblib"

declared = json.loads(config.read_text())
needed = {
    value
    for key, value in declared.items()
    if key.endswith(("_oof", "_test", "_bank", "_definition", "_selections"))
    and isinstance(value, str)
}
needed.update({declared["compounds"], declared["test_mapping"], *declared["complement_models"]})
if needed != set(inputs):
    raise ValueError(f"Input mismatch: missing={needed-set(inputs)}, extra={set(inputs)-needed}")
for path in inputs.values():
    if not path.is_file():
        raise FileNotFoundError(path)

command = [
    "pixi", "run", "cyp", "run",
    "--track", "classification",
    "--purpose", "Fit unrestricted AP Caruana on the weight-guided 90-percent CYP2D6 pool",
    "--config", str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_ap_mass90_fullfit.py"],
    cwd=root,
    check=True,
)
