"""Select the official regression submission columns without changing predictions."""

import os
from pathlib import Path

import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
frame = pl.read_parquet(next(inp.glob("*-tabicl-regression.parquet")))
columns = ["SMILES", "Molecule_Name"] + [
    f"{enzyme}_pIC50_direct_inhibition" for enzyme in ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
]
frame.select(columns).write_parquet(out / "regression.parquet")
