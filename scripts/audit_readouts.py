"""Audit full OOF coverage and score readouts against official ST-RAE."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from split_baselines import metrics

from cyp_challenge.storage import write_json


def main():
    inputs, output = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    spec = importlib.util.spec_from_file_location(
        "official", next(inputs.glob("*-custom_scoring_functions.py"))
    )
    official = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(official)
    predictions = pl.concat([pl.read_parquet(p) for p in sorted(inputs.glob("*-oof.parquet"))])
    summaries = [json.loads(p.read_text()) for p in sorted(inputs.glob("*-summary.json"))]
    assert len(summaries) == 5 and {s["config"]["fold"] for s in summaries} == set(range(5))
    labels = pl.read_parquet(next(inputs.glob("*-measurements.parquet"))).filter(
        pl.col("endpoint").is_in(config["endpoints"])
    )
    assert (
        predictions.select("compound_id", "endpoint").n_unique()
        == predictions.height
        == labels.height
    )
    joined = labels.join(predictions, on=["compound_id", "endpoint"], how="left", validate="1:1")
    assert joined["prediction"].null_count() == 0
    assert (joined["value"] == joined["y_true"]).all()
    scores = []
    for endpoint in config["endpoints"]:
        frame = joined.filter(pl.col("endpoint") == endpoint)
        classification = endpoint.endswith("is_TDI")
        for fold in [-1, 0, 1, 2, 3, 4]:
            f = frame if fold == -1 else frame.filter(pl.col("fold") == fold)
            y, p = f["y_true"].to_numpy(), f["prediction"].to_numpy()
            lower, upper = f["conf_low"].to_numpy(), f["conf_high"].to_numpy()
            m = metrics(y, p, lower, upper, classification)
            if not classification:
                expected = official.rae_soft_threshold_absolute_error(
                    y, p, y_true_upper=upper, y_true_lower=lower
                )
                assert np.isclose(expected, m["ST-RAE"], rtol=1e-12, atol=1e-12)
            if fold >= 0:
                saved = next(s for s in summaries if s["config"]["fold"] == fold)
                row = next(r for r in saved["metrics"] if r["endpoint"] == endpoint)
                assert all(np.isclose(v, row[k], rtol=1e-12, atol=1e-12) for k, v in m.items())
            scores.append(
                {
                    "endpoint": endpoint,
                    "fold": fold,
                    "n": len(y),
                    **m,
                    "predicted_positive_rate": float((p >= 0.5).mean()) if classification else None,
                    "true_positive_rate": float(y.mean()) if classification else None,
                }
            )
    table = pl.DataFrame(scores, infer_schema_length=None)
    table.write_parquet(output / "metrics.parquet")
    table.filter(pl.col("fold") == -1).write_csv(output / "pooled_metrics.csv")
    predictions.write_parquet(output / "oof.parquet")
    macro = (
        table.filter(pl.col("fold") == -1)
        .select(
            pl.col("ST-RAE").mean(),
            pl.col("MAE").mean(),
            pl.col("MCC").mean(),
            pl.col("AUROC").mean(),
        )
        .to_dicts()[0]
    )
    write_json(
        output / "audit.json",
        {
            "passed": True,
            "rows": predictions.height,
            "macro": macro,
            "pooled": table.filter(pl.col("fold") == -1).to_dicts(),
        },
    )
    print(json.dumps(macro), flush=True)
    print(table.filter(pl.col("fold") == -1), flush=True)


if __name__ == "__main__":
    main()
