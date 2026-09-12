"""Bridge author blinded molecules to exact canonical prediction IDs, no blinded labels."""

import json
import os
from pathlib import Path

import polars as pl
from rdkit import Chem

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def path(i):
    return next(inp.glob(f"{i:03d}-*"))


blind = pl.read_csv(path(0), infer_schema_length=None)
pred = pl.read_parquet(path(1))
aliases = pl.read_parquet(path(2))
testmap = pl.read_parquet(path(3))
assert blind.height == 750 and set(blind.columns) == {"molecule_name", "smiles"}
names = aliases.select("molecule_name", "compound_id").unique()
assert names["molecule_name"].n_unique() == names.height
b = blind.join(names, on="molecule_name", validate="1:1")
assert b["compound_id"].null_count() == 0
assert set(b["compound_id"]) == set(testmap["compound_id"])
comparison = b.join(
    pred.select("compound_id", pl.col("smiles").alias("previous_smiles")),
    on="compound_id",
    validate="1:1",
)
changed = comparison.filter(pl.col("smiles") != pl.col("previous_smiles"))
changed.write_csv(out / "changed-prediction-smiles.csv")
new = (
    pred.join(
        b.select("compound_id", pl.col("smiles").alias("blind_smiles")),
        on="compound_id",
        how="left",
        validate="1:1",
    )
    .with_columns(pl.coalesce("blind_smiles", "smiles").alias("smiles"))
    .drop("blind_smiles")
    .with_columns(
        pl.col("smiles")
        .map_elements(
            lambda s: Chem.MolToInchiKey(Chem.MolFromSmiles(s)).split("-")[0],
            return_dtype=pl.String,
        )
        .alias("connectivity")
    )
)
assert new.height == 6896 and new["compound_id"].n_unique() == 6896
new.write_parquet(out / "prediction-compounds.parquet")
(out / "verification.json").write_text(
    json.dumps(
        dict(
            blinded750=True,
            no_blinded_labels=True,
            prediction_smiles_changes=changed.height,
            canonical_ids_exact=True,
        ),
        indent=2,
    )
)
print("changed", changed.height)
