"""Package the retained submitted D6 vector with the audited A4 TabICL Caruana result."""

import json
import os
from pathlib import Path

import polars as pl

from cyp_challenge.tracks import ENDPOINTS, validate

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
reference = pl.read_csv(inp / cfg["reference"])
aliases = (
    pl.read_parquet(inp / cfg["aliases"])
    .filter(pl.col("role") == "test")
    .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
    .unique()
)
identities = reference.select("SMILES", "Molecule_Name").join(
    aliases, on="Molecule_Name", validate="1:1"
)
d6 = (
    pl.read_parquet(inp / cfg["current"])
    .filter(
        (pl.col("endpoint") == "CYP2D6_is_TDI")
        & (pl.col("method") == "all_plain_pfn_graph_mean")
    )
    .select("compound_id", pl.col("binary").alias("CYP2D6_is_TDI"))
)
a4 = (
    pl.read_parquet(inp / cfg["caruana"])
    .filter(
        (pl.col("endpoint") == "CYP3A4_is_TDI")
        & (pl.col("method") == "caruana_auc")
        & (pl.col("fold") == -1)
    )
    .select("compound_id", pl.col("binary").alias("CYP3A4_is_TDI"))
)
candidate = (
    identities.join(d6, on="compound_id", validate="1:1")
    .join(a4, on="compound_id", validate="1:1")
    .select("SMILES", "Molecule_Name", *ENDPOINTS["classification"])
    .sort("Molecule_Name")
)
assert candidate.height == 750
candidate.write_parquet(out / "retain-d6-tabicl-caruana-a4.parquet")
candidate.write_csv(out / "retain-d6-tabicl-caruana-a4.csv")
report = dict(
    validation=validate(
        out / "retain-d6-tabicl-caruana-a4.parquet", "classification", inp / cfg["reference"]
    ),
    positives={ep: int(candidate[ep].sum()) for ep in ENDPOINTS["classification"]},
    source_methods={
        "CYP2D6_is_TDI": "all_plain_pfn_graph_mean",
        "CYP3A4_is_TDI": "caruana_auc",
    },
)
(out / "summary.json").write_text(json.dumps(report, indent=2))
print(json.dumps(report))
