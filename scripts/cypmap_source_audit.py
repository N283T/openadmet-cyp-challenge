"""Audit CYP-MAP source molecules against challenge connectivity."""

import json
import os
from pathlib import Path

import polars as pl
from rdkit import Chem

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def connectivity(molecule: Chem.Mol | None) -> str | None:
    if molecule is None:
        return None
    key = Chem.MolToInchiKey(molecule)
    return key.split("-", 1)[0] if key else None


compounds = pl.read_parquet(inp / cfg["compounds"])
test_ids = set(pl.read_parquet(inp / cfg["test_mapping"])["compound_id"])
challenge = pl.DataFrame(
    {
        "compound_id": compound_id,
        "connectivity": connectivity(Chem.MolFromSmiles(smiles)),
        "challenge_role": "test" if compound_id in test_ids else "train",
    }
    for compound_id, smiles in compounds.select("compound_id", "standardized_smiles").iter_rows()
)
challenge = challenge.group_by("connectivity").agg(
    pl.col("compound_id").sort().str.join("|").alias("compound_ids"),
    pl.col("challenge_role").eq("train").any().alias("overlaps_train"),
    pl.col("challenge_role").eq("test").any().alias("overlaps_test"),
)

records = []
invalid = []
for source_split, filename in cfg["sdf_files"].items():
    supplier = Chem.SDMolSupplier(str(inp / filename), removeHs=False)
    for source_row, molecule in enumerate(supplier):
        if molecule is None:
            invalid.append({"source_split": source_split, "source_row": source_row})
            continue
        records.append(
            {
                "source_split": source_split,
                "source_row": source_row,
                "source_id": molecule.GetProp("ID") if molecule.HasProp("ID") else None,
                "source_name": molecule.GetProp("Name") if molecule.HasProp("Name") else None,
                "connectivity": connectivity(molecule),
            }
        )

external = pl.DataFrame(records).join(challenge, on="connectivity", how="left", validate="m:1")
overlap = external.filter(pl.col("compound_ids").is_not_null()).sort("source_split", "compound_ids")
external.write_parquet(out / "cypmap-source-records.parquet")
overlap.write_parquet(out / "cypmap-challenge-overlap.parquet")
pl.DataFrame(
    invalid,
    schema={"source_split": pl.String, "source_row": pl.Int64},
).write_parquet(out / "cypmap-invalid-records.parquet")

summary = {
    "source_commit": cfg["source_commit"],
    "instances": external.height,
    "unique_connectivities": external["connectivity"].n_unique(),
    "invalid_records": len(invalid),
    "challenge_overlap_instances": overlap.height,
    "challenge_overlap_connectivities": overlap["connectivity"].n_unique(),
    "challenge_train_overlap_connectivities": overlap.filter(
        pl.col("overlaps_train").fill_null(False)
    )["connectivity"].n_unique(),
    "challenge_test_overlap_connectivities": overlap.filter(
        pl.col("overlaps_test").fill_null(False)
    )["connectivity"].n_unique(),
    "interpretation": (
        "Frozen published checkpoints may be evaluated only as external features; "
        "overlap counts qualify validation and prohibit treating the readout as a clean "
        "external benchmark."
    ),
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
