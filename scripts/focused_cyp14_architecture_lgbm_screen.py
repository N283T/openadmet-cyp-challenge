"""Run matched raw LightGBM readouts for already-trained focused CYP encoders."""

import importlib.util
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import polars as pl

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
root = Path(__file__).resolve().parents[1]

spec = importlib.util.spec_from_file_location("official", inp / "custom_scoring_functions.py")
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)

rows = []
for variant in cfg["variants"]:
    variant_out = out / variant["name"]
    variant_out.mkdir()
    variant_config = {
        "method": variant["name"],
        "backend": "lgbm",
        "feature_file": variant["file"],
        "feature_prefix": cfg["feature_prefix"],
        "pca": 0,
        "mask_observed_validation": False,
        "split_run": cfg["split_run"],
        "endpoints": cfg["endpoints"],
    }
    variant_config_path = variant_out / "config.json"
    variant_config_path.write_text(json.dumps(variant_config, indent=2) + "\n")
    env = {
        **os.environ,
        "CYP_OUTPUT_DIR": str(variant_out),
        "CYP_CONFIG_PATH": str(variant_config_path),
    }
    subprocess.run(
        ["python", "scripts/external_embedding_readout.py"], cwd=root, env=env, check=True
    )

all_variants = cfg["variants"] + cfg["controls"]
for variant in all_variants:
    path = (
        out / variant["name"] / "external-oof.parquet"
        if variant in cfg["variants"]
        else inp / variant["file"]
    )
    frame = pl.read_parquet(path)
    for endpoint in cfg["endpoints"]:
        selected = frame.filter(pl.col("endpoint") == endpoint)
        value = float(
            scorer.rae_soft_threshold_absolute_error(
                selected["y_true"].to_numpy(),
                selected["prediction"].to_numpy(),
                y_true_lower=selected["conf_low"].to_numpy(),
                y_true_upper=selected["conf_high"].to_numpy(),
            )
        )
        rows.append({"method": variant["name"], "endpoint": endpoint, "ST_RAE": value})

metrics = pl.DataFrame(rows).sort("endpoint", "ST_RAE")
macro = metrics.group_by("method").agg(pl.col("ST_RAE").mean()).sort("ST_RAE")
assert np.isfinite(metrics["ST_RAE"].to_numpy()).all()
metrics.write_csv(out / "metrics.csv")
macro.write_csv(out / "macro.csv")
(out / "summary.json").write_text(
    json.dumps(
        {"config": cfg, "metrics": metrics.to_dicts(), "macro": macro.to_dicts()}, indent=2
    )
    + "\n"
)
print(metrics, flush=True)
print(macro, flush=True)
