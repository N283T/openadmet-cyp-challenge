"""Independently verify target-wise component pruning and packaged predictions."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl

from cyp_challenge.tracks import validate

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
keys = ["compound_id", "endpoint", "fold"]

spec = importlib.util.spec_from_file_location("official_scorer", inp / cfg["scorer"])
assert spec is not None and spec.loader is not None
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)


def score(frame: pl.DataFrame, prediction: np.ndarray) -> float:
    return float(
        scorer.rae_soft_threshold_absolute_error(
            frame["y_true"].to_numpy(),
            prediction,
            y_true_lower=frame["conf_low"].to_numpy(),
            y_true_upper=frame["conf_high"].to_numpy(),
        )
    )


crossfit = pl.read_csv(inp / cfg["crossfit_scores"])
full_best = pl.read_csv(inp / cfg["full_oof_best"])
oof = {}
test = {}
for component in cfg["components"]:
    name = component["name"]
    oof[name] = pl.read_parquet(inp / component["oof"]).sort(keys)
    test[name] = pl.read_parquet(inp / component["test"]).sort(keys)

base = next(iter(oof.values()))
for frame in oof.values():
    for column in keys + ["y_true", "conf_low", "conf_high", "is_neighbor"]:
        assert frame[column].equals(base[column]), (column, frame)

packaged = pl.read_parquet(inp / cfg["packaged_long"]).sort(keys)
assert packaged.height == 3000 and not packaged.select(keys).is_duplicated().any()
assert validate(inp / cfg["packaged_wide"], "regression", inp / cfg["reference"])["valid"]

rows = []
for endpoint, selected in cfg["endpoint_components"].items():
    count = len(selected)
    endpoint_crossfit = crossfit.filter(pl.col("endpoint") == endpoint).sort("component_count")
    eligible = endpoint_crossfit.filter(
        (pl.col("delta_ci_low") <= 0) & (pl.col("delta_ci_high") >= 0)
    )
    assert eligible["component_count"].min() == count
    fitted = full_best.filter(
        (pl.col("endpoint") == endpoint) & (pl.col("component_count") == count)
    ).row(0, named=True)
    assert fitted["components"].split("+") == selected

    endpoint_oof = base.filter(pl.col("endpoint") == endpoint)
    p_oof = np.mean(
        np.stack(
            [oof[name].filter(pl.col("endpoint") == endpoint)["prediction"] for name in selected]
        ),
        axis=0,
    )
    oof_score = score(endpoint_oof, p_oof)
    assert abs(oof_score - fitted["ST_RAE"]) < 1e-12

    endpoint_test = packaged.filter(pl.col("endpoint") == endpoint)
    p_test = np.mean(
        np.stack(
            [test[name].filter(pl.col("endpoint") == endpoint)["prediction"] for name in selected]
        ),
        axis=0,
    )
    assert np.allclose(endpoint_test["prediction"], p_test, rtol=0, atol=1e-12)
    rows.append(
        {
            "endpoint": endpoint,
            "component_count": count,
            "crossfit_ST_RAE": endpoint_crossfit.filter(pl.col("component_count") == count)[
                "ST_RAE"
            ].item(),
            "full_oof_ST_RAE": oof_score,
            "max_test_reconstruction_error": float(
                np.max(np.abs(endpoint_test["prediction"].to_numpy() - p_test))
            ),
        }
    )

scores = pl.DataFrame(rows).sort("endpoint")
scores.write_csv(out / "scores.csv")
verification = {
    "status": "passed",
    "endpoint_component_count": int(scores["component_count"].sum()),
    "previous_endpoint_component_count": 20,
    "crossfit_macro_ST_RAE": float(scores["crossfit_ST_RAE"].mean()),
    "full_oof_macro_ST_RAE": float(scores["full_oof_ST_RAE"].mean()),
    "packaged_rows": packaged.height,
    "selection_rule_verified": True,
    "test_means_reconstructed": True,
}
(out / "verification.json").write_text(json.dumps(verification, indent=2))
print(json.dumps(verification, indent=2), flush=True)
