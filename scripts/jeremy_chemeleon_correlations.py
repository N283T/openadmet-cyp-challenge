"""Matched OOF correlations and fixed equal-blend metrics; no weight fitting."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
keys = ["compound_id", "endpoint", "fold"]
base = pl.read_parquet(inp / cfg["new_oof"]).sort(keys)
assert base.height == 6384 and base.select(keys).n_unique() == 6384
rows = []
for method, filename in cfg["controls"].items():
    other = pl.read_parquet(inp / filename).sort(keys)
    assert other.height == 6384 and other.select(keys).n_unique() == 6384
    for column in keys + ["y_true", "is_neighbor", "conf_low", "conf_high"]:
        assert base[column].equals(other[column]), column
    joined = base.with_columns(other["prediction"].alias("control_prediction"))
    for (endpoint,), group in joined.group_by("endpoint"):
        for subset in ["all", "neighbors"]:
            data = group if subset == "all" else group.filter(pl.col("is_neighbor"))
            y = data["y_true"].to_numpy()
            new, control = data["prediction"].to_numpy(), data["control_prediction"].to_numpy()
            rows.append(
                dict(
                    endpoint=endpoint,
                    subset=subset,
                    control=method,
                    n=data.height,
                    prediction_pearson=float(np.corrcoef(new, control)[0, 1]),
                    residual_pearson=float(np.corrcoef(new - y, control - y)[0, 1]),
                )
            )
correlations = pl.DataFrame(rows).sort("control", "subset", "endpoint")
correlations.write_csv(out / "correlations.csv")
scores = pl.concat(
    [
        pl.read_csv(inp / item["file"]).filter(pl.col("method") == item["method"])
        for item in cfg["scores"]
    ]
)
scores.write_csv(out / "scores.csv")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            config=cfg,
            scores=scores.to_dicts(),
            correlations=correlations.to_dicts(),
            weights_fitted=False,
            matched_oof=base.height,
        ),
        indent=2,
    )
)
print(scores)
print(correlations)
