"""Extract human NVS CYP assay observations from official EPA release."""

import json
import os
from pathlib import Path

import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
t = pl.read_csv(
    next(inp.glob("*-EXPORT_LVL5&6_ASID5_NVS_200730.csv")), null_values="NA", infer_schema=False
)
nvs = t.filter(pl.col("aenm").str.contains(r"^NVS_ADME_hCYP(1A2|2C19|2C9|2D6|3A4)$"))
nvs.write_parquet(out / "nvs-raw.parquet")
nvs.select("chid", "casn", "chnm", "code").unique().write_csv(out / "nvs-identities.csv")
counts = nvs.group_by("aeid", "aenm", "hitc", "tested_conc_unit").agg(
    pl.len().alias("n"), pl.col("code").n_unique().alias("codes")
)
counts.write_csv(out / "counts.csv")
print(nvs.select("chid", "casn", "chnm", "code").unique().head(15), flush=True)
print(counts, flush=True)
(out / "summary.json").write_text(
    json.dumps(
        {"rows": nvs.height, "columns": nvs.columns, "identities": nvs["code"].n_unique()}, indent=2
    )
)
