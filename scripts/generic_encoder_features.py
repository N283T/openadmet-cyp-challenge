"""Extract one generic frozen encoder feature table through its isolated runtime."""

import json
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import polars as pl

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
scratch = Path(os.environ["TMPDIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def input_at(index: int) -> Path:
    paths = list(inp.glob(f"{index:03d}-*"))
    if len(paths) != 1:
        raise ValueError(f"Expected one input at {index}, found {paths}")
    return paths[0]


compounds = pl.read_parquet(input_at(0)).sort("compound_id")
if cfg["limit"]:
    compounds = compounds.head(cfg["limit"])
csv_path = scratch / "smiles.csv"
compounds.select("standardized_smiles").rename({"standardized_smiles": "smiles"}).write_csv(
    csv_path
)
npz_path = scratch / f"{cfg['encoder']}.npz"

if cfg["encoder"] == "kermt_v2":
    checkpoint = input_at(1)
    command = [
        "pixi",
        "run",
        "--manifest-path",
        cfg["runtime_manifest"],
        "python",
        cfg["driver"],
        "fingerprint",
        "--data_path",
        str(csv_path),
        "--checkpoint_path",
        str(checkpoint),
        "--output_path",
        str(npz_path),
        "--fingerprint_source",
        "both",
        "--batch_size",
        str(cfg["batch_size"]),
    ]
    environment = os.environ.copy()
    environment["PYTHONPATH"] = cfg["source_dir"]
    environment["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
elif cfg["encoder"] == "molformer_c3":
    model_dir = scratch / "molformer-model"
    model_dir.mkdir()
    for index, name in enumerate(cfg["model_files"], start=1):
        shutil.copy2(input_at(index), model_dir / name)
    code_start = 1 + len(cfg["model_files"])
    for offset, name in enumerate(cfg["code_files"]):
        shutil.copy2(input_at(code_start + offset), model_dir / name)
    for metadata_name in ("config.json", "tokenizer_config.json"):
        metadata_path = model_dir / metadata_name
        metadata_path.chmod(0o600)
        metadata = json.loads(metadata_path.read_text())
        auto_map = metadata.get("auto_map", {})
        for key, value in auto_map.items():
            if isinstance(value, str):
                auto_map[key] = value.split("--", maxsplit=1)[-1]
            else:
                auto_map[key] = [item.split("--", maxsplit=1)[-1] for item in value]
        metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    command = [
        "pixi",
        "run",
        "--manifest-path",
        cfg["runtime_manifest"],
        "python",
        cfg["driver"],
        "--data-path",
        str(csv_path),
        "--model-path",
        str(model_dir),
        "--output-path",
        str(npz_path),
        "--batch-size",
        str(cfg["batch_size"]),
        "--pooling",
        cfg["pooling_mode"],
    ]
    environment = os.environ.copy()
else:
    raise ValueError(cfg["encoder"])

subprocess.run(command, check=True, env=environment)
with np.load(npz_path) as archive:
    features = archive["fps"].astype(np.float32)
if features.shape[0] != compounds.height:
    raise ValueError(f"Feature row mismatch: {features.shape[0]} != {compounds.height}")
prefix = cfg["feature_prefix"]
feature_frame = pl.DataFrame(
    features, schema=[f"{prefix}_{index:04d}" for index in range(features.shape[1])]
).insert_column(0, compounds["compound_id"])
feature_frame.write_parquet(out / "features.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "config": cfg,
            "rows": feature_frame.height,
            "dimensions": features.shape[1],
            "finite": bool(np.isfinite(features).all()),
            "mean_absolute": float(np.abs(features).mean()),
            "mean_column_std": float(features.std(axis=0).mean()),
        },
        indent=2,
    )
    + "\n"
)
print(json.dumps({"rows": feature_frame.height, "dimensions": features.shape[1]}), flush=True)
