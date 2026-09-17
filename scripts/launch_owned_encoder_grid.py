"""Train project-owned external encoders and make CYP2D6 atomic OOF readouts."""

import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GRID = json.loads((ROOT / "configs/external-training/owned-encoder-grid-v1.json").read_text())
DATASET = ROOT / f"artifacts/runs/shared/{GRID['dataset_run']}/outputs"
CHALLENGE = ROOT / "artifacts/runs/shared/20260909T132734119233Z-a7063ccd/outputs"
ASSIGNMENTS = (
    ROOT
    / "artifacts/runs/regression/20260911T104703655337Z-a3923255/outputs"
    / "parent-assignments.parquet"
)
TABPFN = Path.home() / ".cache/tabpfn/tabpfn-v2.6-regressor-v2.6_default.ckpt"
OUTPUT = ROOT / "artifacts/launches/owned_encoder_grid_v1_manifest.json"


def run_recorded(
    *, track: str, purpose: str, config: Path, inputs: dict[str, Path], script: str
) -> str:
    command = [
        "pixi",
        "run",
        "cyp",
        "run",
        "--track",
        track,
        "--purpose",
        purpose,
        "--config",
        str(config),
    ]
    for alias, path in inputs.items():
        command.extend(["--input", f"{alias}={path}"])
    result = subprocess.run(
        command + ["--", "python", script], cwd=ROOT, text=True, capture_output=True
    )
    print(result.stdout, flush=True)
    if result.stderr:
        print(result.stderr, flush=True)
    result.check_returncode()
    run_ids = re.findall(r"run_id=(\S+)", result.stdout)
    if len(run_ids) != 1:
        raise ValueError(f"Expected one run ID; got {run_ids}")
    return run_ids[0]


config_dir = ROOT / "configs/owned-encoder-grid"
config_dir.mkdir(exist_ok=True)
manifest = []
for width in GRID["widths"]:
    for seed in GRID["seeds"]:
        key = f"{width}:{seed}"
        encoder_run = GRID.get("reuse", {}).get(key)
        if encoder_run is None:
            tag = f"owned_external_w{width}_s{seed}_v1"
            encoder_config = {
                key: value
                for key, value in GRID.items()
                if key
                not in {
                    "widths",
                    "seeds",
                    "reuse",
                    "readouts",
                    "pca_dimensions",
                    "endpoints",
                    "protocol",
                }
            }
            encoder_config.update(
                {
                    "architecture": "owned_width_grid",
                    "width": width,
                    "seed": seed,
                    "purpose": GRID["protocol"],
                }
            )
            encoder_config_path = config_dir / f"{tag}.json"
            encoder_config_path.write_text(json.dumps(encoder_config, indent=2) + "\n")
            encoder_run = run_recorded(
                track="shared",
                purpose=f"Project-owned external encoder width {width}, training seed {seed}",
                config=encoder_config_path,
                inputs={
                    "external-compounds.parquet": DATASET / "external-compounds.parquet",
                    "observations.parquet": DATASET / "observations.parquet",
                    "compounds.parquet": CHALLENGE / "compounds.parquet",
                },
                script="scripts/train_external_medium.py",
            )
        feature_path = (
            ROOT / f"artifacts/runs/shared/{encoder_run}/outputs/external-embeddings.parquet"
        )
        for backend in GRID["readouts"]:
            for pca in GRID["pca_dimensions"]:
                representation = f"pca{pca}" if pca else f"raw{width}"
                tag = f"owned_external_w{width}_s{seed}_{representation}_{backend}_d6_v1"
                readout_config = {
                    "method": tag,
                    "backend": backend,
                    "feature_file": "external-embeddings.parquet",
                    "feature_prefix": "emb_",
                    "pca": pca,
                    "mask_observed_validation": False,
                    "split_run": "20260911T104703655337Z-a3923255",
                    "endpoints": GRID["endpoints"],
                    "encoder_run": encoder_run,
                    "encoder_width": width,
                    "encoder_seed": seed,
                }
                inputs = {
                    "external-embeddings.parquet": feature_path,
                    "measurements.parquet": CHALLENGE / "measurements.parquet",
                    "parent-assignments.parquet": ASSIGNMENTS,
                }
                if backend == "tabpfn":
                    readout_config["weights"] = "tabpfn-v2.6-regressor-v2.6_default.ckpt"
                    inputs["tabpfn-v2.6-regressor-v2.6_default.ckpt"] = TABPFN
                readout_config_path = config_dir / f"{tag}.json"
                readout_config_path.write_text(json.dumps(readout_config, indent=2) + "\n")
                readout_run = run_recorded(
                    track="regression",
                    purpose=(
                        f"CYP2D6 atomic {backend} readout of owned width {width} "
                        f"seed {seed}, {representation}"
                    ),
                    config=readout_config_path,
                    inputs=inputs,
                    script="scripts/external_embedding_readout.py",
                )
                manifest.append(
                    {
                        "name": tag,
                        "width": width,
                        "seed": seed,
                        "pca": pca,
                        "backend": backend,
                        "encoder_run": encoder_run,
                        "readout_run": readout_run,
                        "oof": str(
                            ROOT
                            / f"artifacts/runs/regression/{readout_run}/outputs"
                            / "external-oof.parquet"
                        ),
                    }
                )
                OUTPUT.write_text(json.dumps(manifest, indent=2) + "\n")

OUTPUT.write_text(json.dumps(manifest, indent=2) + "\n")
print(json.dumps({"members": len(manifest), "manifest": manifest}, indent=2), flush=True)
