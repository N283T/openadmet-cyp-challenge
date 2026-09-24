"""Launch blind decision stability for the five AP-mass90 outer ensembles."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs/classification"
fixed = runs / "20260913T123045515718Z-0df22d54"
config = root / "configs/analysis/d6-ap-mass90-blind-stability-v1.json"
inputs = {
    "outer-details.parquet": runs / "20260923T171350543926Z-9081bf31/outputs/outer-details.parquet",
    "decision-test.parquet": runs / "20260923T173655387004Z-97d71dce/outputs/decision-test.parquet",
    "support-test.parquet": runs / "20260918T085712373599Z-8829d528/outputs/test.parquet",
    "lgbm-test.parquet": runs / "20260919T155849130417Z-0394b600/outputs/test.parquet",
    "readout-test.parquet": runs / "20260923T161449619218Z-a5f9f7ea/outputs/test.parquet",
    "tabpfn-test.parquet": runs / "20260923T161718040624Z-8a9676e7/outputs/test.parquet",
    "mace-test.parquet": runs / "20260914T051014577658Z-7b1b64c2/outputs/test.parquet",
    "activity-test.parquet": runs / "20260915T065610180920Z-f532943e/outputs/test.parquet",
    "mechanism-test.parquet": runs / "20260914T063614433212Z-c0dc3f32/outputs/test.parquet",
    "pfn-pca-test.parquet": runs
    / "20260912T121146867140Z-b0a05986/inputs/005-test-fold-probabilities.parquet",
    "pfn-raw-test.parquet": runs
    / "20260915T112944344816Z-69abecf5/outputs/test-probabilities.parquet",
    "pool-definition.json": root / "configs/classification/d6-flat-common-pool-caruana-v1.json",
    "compounds.parquet": root
    / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
    "test-mapping.parquet": fixed / "inputs/002-test_mapping.parquet",
}
for fold in range(5):
    inputs[f"complement-fold{fold}.joblib"] = fixed / f"outputs/complement_nb_a1-fold{fold}.joblib"
declared = json.loads(config.read_text())
needed = {
    declared[key]
    for key in [
        "outer_details", "decision_test", "support_test", "lgbm_test",
        "readout_test", "tabpfn_test", "mace_test", "activity_test",
        "mechanism_test", "pfn_pca_test", "pfn_raw_test", "pool_definition",
        "compounds", "test_mapping",
    ]
}
needed.update(declared["complement_models"])
if needed != set(inputs):
    raise ValueError(f"Input mismatch: missing={needed-set(inputs)}, extra={set(inputs)-needed}")
for path in inputs.values():
    if not path.is_file():
        raise FileNotFoundError(path)
command = [
    "pixi", "run", "cyp", "run", "--track", "classification",
    "--purpose", "Audit blind call stability across five AP-mass90 outer ensembles",
    "--config", str(config),
]
for alias, path in inputs.items():
    command.extend(["--input", f"{alias}={path}"])
subprocess.run(
    command + ["--", "python", "scripts/classification_d6_ap_mass90_blind_stability.py"],
    cwd=root, check=True,
)
