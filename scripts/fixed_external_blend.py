"""Key-checked fixed equal OOF blend; no fitting or weight search."""

import json
import os
from pathlib import Path

import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
keys = ["compound_id", "endpoint", "fold"]
base, predictions = None, []
for component in cfg["components"]:
    frame = pl.read_parquet(inp / component["file"]).sort(keys)
    assert frame.height == 6384 and not frame.select(keys).is_duplicated().any()
    if base is None:
        base = frame
    for column in keys + ["y_true", "conf_low", "conf_high", "is_neighbor"]:
        assert frame[column].equals(base[column]), column
    predictions.append(frame[component["column"]])
prediction = sum(predictions) / len(predictions)
assert prediction.is_finite().all()
base.select(keys + ["y_true", "conf_low", "conf_high", "is_neighbor"]).with_columns(
    prediction.alias("prediction"),
    pl.lit(cfg["method"]).alias("method"),
    pl.lit(cfg["split_run"]).alias("split_id"),
).write_parquet(out / "external-oof.parquet")
(out / "summary.json").write_text(json.dumps(cfg, indent=2))
