"""Append shifted 3A4 TDI-condition pIC50 of TDI-negative compounds as direct labels.

Full-data counterpart of `regression_3a4_tdi_proxy.py`: the shift is the median
TDI-minus-direct difference over every TDI-negative compound carrying both labels,
and interval bounds move with the value. Original measurements are unchanged.
"""

import json
import os
from pathlib import Path

import polars as pl

DIRECT = "CYP3A4_pIC50_direct_inhibition"
TDI = "CYP3A4_pIC50_TDI_condition"
LABEL = "CYP3A4_is_TDI"


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    measurements = pl.read_parquet(inp / "measurements.parquet")
    wide = measurements.filter(pl.col("endpoint").is_in([DIRECT, TDI, LABEL])).pivot(
        on="endpoint", index="compound_id", values="value", aggregate_function="first"
    )
    negative = pl.col(LABEL) == 0
    paired = wide.filter(pl.col(DIRECT).is_not_null() & pl.col(TDI).is_not_null() & negative)
    shift = float((paired[TDI] - paired[DIRECT]).median())
    proxy_ids = wide.filter(pl.col(DIRECT).is_null() & pl.col(TDI).is_not_null() & negative)[
        "compound_id"
    ]
    assert proxy_ids.len() == cfg["expected_proxy_rows"]
    tdi_rows = measurements.filter(
        (pl.col("endpoint") == TDI) & pl.col("compound_id").is_in(proxy_ids.to_list())
    )
    assert tdi_rows["compound_id"].n_unique() == tdi_rows.height == proxy_ids.len()
    proxy = tdi_rows.select(
        ("proxy-tdi-" + pl.col("compound_id")).alias("measurement_id"),
        "compound_id",
        pl.lit(DIRECT).alias("endpoint"),
        pl.col("value") - shift,
        pl.col("conf_low") - shift,
        pl.col("conf_high") - shift,
        "std",
    )
    augmented = pl.concat([measurements, proxy.cast(measurements.schema)])
    assert (
        augmented.select("compound_id", "endpoint")
        .filter(pl.col("endpoint") == DIRECT)
        .is_unique()
        .all()
    )
    augmented.write_parquet(out / "measurements.parquet")
    proxy.write_parquet(out / "proxy-rows.parquet")
    summary = {
        "shift": shift,
        "paired_tdi_negative": paired.height,
        "proxy_rows": proxy.height,
        "direct_rows_before": measurements.filter(pl.col("endpoint") == DIRECT).height,
        "direct_rows_after": augmented.filter(pl.col("endpoint") == DIRECT).height,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
