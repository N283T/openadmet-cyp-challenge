"""Record matched frozen upstream medium/CheMeleon control scores and provenance."""

import json
import os
from pathlib import Path

import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
tables = []
for representation, item in cfg["representations"].items():
    table = pl.read_csv(inp / item["file"]).filter(pl.col("method") == item["method"])
    assert table.height > 0
    tables.append(table.with_columns(pl.lit(representation).alias("representation")))
scores = pl.concat(tables)
scores.write_csv(out / "matched-frozen-control.csv")
summary = json.loads((inp / cfg["export_summary"]).read_text())
assert summary["dimensions"] == 2048 and summary["pretraining_output_heads"] == 47
assert summary["challenge_labels_used_during_extraction"] is False
result = dict(export=summary, matched_scores=scores.to_dicts())
(out / "summary.json").write_text(json.dumps(result, indent=2))
print(json.dumps(result, indent=2))
