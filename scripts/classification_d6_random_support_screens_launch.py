"""Launch random-fold D6 support-mismatch diagnostics."""

import json
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[1]
runs = root / "artifacts/runs"
views = runs / "regression/20260917T111518065286Z-801ab060/outputs"


def launch(config_name: str, purpose: str, script: str, inputs: dict[str, Path]) -> None:
    config = root / "configs/classification" / config_name
    declared = json.loads(config.read_text())
    if "blocks" in declared:
        needed = {spec["file"] for spec in declared["blocks"].values()}
    else:
        needed = set()
    for key in [
        "compounds",
        "labels",
        "test_mapping",
        "assignments",
        "current_bank",
        "tabpfn_weights",
    ]:
        if key in declared:
            needed.add(declared[key])
    if needed != set(inputs):
        raise ValueError(
            f"Input mismatch for {config_name}: missing={needed - set(inputs)}, "
            f"extra={set(inputs) - needed}"
        )
    for path in inputs.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        "classification",
        "--purpose",
        purpose,
        "--config",
        str(config),
    ]
    for alias, path in inputs.items():
        command.extend(["--input", f"{alias}={path}"])
    subprocess.run(command + ["--", "python", script], cwd=root, check=True)


shared = {
    "test-mapping.parquet": runs
    / "shared/20260909T132734119233Z-a7063ccd/outputs/test_mapping.parquet",
    "assignments.parquet": runs
    / "shared/20260909T135130268493Z-edba6fc0/outputs/assignments.parquet",
    "current-prediction-bank.parquet": runs
    / "classification/20260915T121009912691Z-73198c15/outputs/prediction-bank.parquet",
}
launch(
    "d6-tabpfn35-fast-random-support-v1.json",
    "Test TabPFN 3.5-fast CYP2D6 under the more test-like random support split",
    "scripts/classification_d6_tabpfn35_fast_feature_screen.py",
    shared
    | {
        "measurements.parquet": runs
        / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
        "tabpfn-v3.5-fast-20260909.safetensors": Path.home()
        / ".cache/tabpfn/tabpfn-v3.5-fast-20260909.safetensors",
        "fp-mordred.parquet": runs
        / "shared/20260912T061043373484Z-2e344d7a/outputs/fp-mordred-features.parquet",
        "jeremy-raw.parquet": runs
        / "shared/20260912T055819206227Z-c060d197/outputs/external-embeddings.parquet",
        "jeremy-scores.parquet": runs
        / "shared/20260912T072105522912Z-53eb675c/outputs/external-scores.parquet",
        "anchor-raw-scores.parquet": views / "jeremy_init_task_anchor-raw-scores.parquet",
        "pub030-raw.parquet": views / "jeremy_pubchem_anchor030-raw.parquet",
        "owned-raw-scores.parquet": views / "owned_init_hybrid-raw-scores.parquet",
        "scratch-raw.parquet": views / "scratch_teacher_only-raw.parquet",
    },
)
launch(
    "d6-vaishnavi-random-support-v1.json",
    "Test the public shallow tree CYP2D6 recipe under the more test-like random split",
    "scripts/classification_d6_vaishnavi_tree_ensemble.py",
    shared
    | {
        "compounds.parquet": runs
        / "shared/20260909T132734119233Z-a7063ccd/outputs/compounds.parquet",
        "classification-labels.parquet": runs
        / "shared/20260909T132734119233Z-a7063ccd/outputs/measurements.parquet",
    },
)
