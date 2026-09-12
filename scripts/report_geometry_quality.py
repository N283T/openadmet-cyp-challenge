"""Audit geometry missingness and standardized molecule connectivity."""

import json
import os
from pathlib import Path

import polars as pl
from rdkit import Chem

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
r, f, c = [pl.read_parquet(next(inp.glob(f"{i:03d}-*"))) for i in range(3)]
counts = r.group_by("compound_id").agg((pl.col("status") == "ok").sum().alias("successful_seeds"))
fragments = pl.DataFrame(
    [
        dict(
            compound_id=x["compound_id"],
            n_fragments=len(Chem.GetMolFrags(Chem.MolFromSmiles(x["standardized_smiles"]))),
        )
        for x in c.iter_rows(named=True)
    ]
)
fragments.write_parquet(out / "fragment-counts.parquet")
quality = dict(
    success_distribution=counts.group_by("successful_seeds")
    .len()
    .sort("successful_seeds")
    .to_dicts(),
    disconnected=fragments.filter(pl.col("n_fragments") > 1).to_dicts(),
    missing=[
        dict(
            feature=k,
            missing=int(f[k].is_null().sum() + f[k].is_nan().sum()),
            nonfinite=int((~f[k].is_finite()).sum()),
        )
        for k in f.columns
        if k != "compound_id"
    ],
)
(out / "quality.json").write_text(json.dumps(quality, indent=2))
print(json.dumps(quality), flush=True)
