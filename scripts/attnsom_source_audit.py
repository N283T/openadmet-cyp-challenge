"""Audit ATTNSOM's atom-level external CYP data against challenge connectivity."""

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
test_mapping = pl.read_parquet(inp / cfg["test_mapping"])
test_ids = set(test_mapping["compound_id"])
challenge_rows = []
for compound_id, smiles in compounds.iter_rows():
    challenge_rows.append(
        {
            "compound_id": compound_id,
            "connectivity": connectivity(Chem.MolFromSmiles(smiles)),
            "challenge_role": "test" if compound_id in test_ids else "train",
        }
    )
challenge = pl.DataFrame(challenge_rows)
if challenge["connectivity"].null_count():
    raise ValueError("Challenge connectivity keys are incomplete")
challenge = challenge.group_by("connectivity").agg(
    pl.col("compound_id").sort().str.join("|").alias("compound_ids"),
    pl.col("challenge_role").eq("train").any().alias("overlaps_train"),
    pl.col("challenge_role").eq("test").any().alias("overlaps_test"),
).with_columns(
    pl.when(pl.col("overlaps_train") & pl.col("overlaps_test"))
    .then(pl.lit("both"))
    .when(pl.col("overlaps_test"))
    .then(pl.lit("test"))
    .otherwise(pl.lit("train"))
    .alias("challenge_role")
)

rows = []
invalid = []
for isoform, filename in cfg["sdf_files"].items():
    supplier = Chem.SDMolSupplier(str(inp / filename), removeHs=False)
    for source_row, molecule in enumerate(supplier):
        if molecule is None:
            invalid.append({"isoform": isoform, "source_row": source_row})
            continue
        primary_som = molecule.GetProp("PRIMARY_SOM") if molecule.HasProp("PRIMARY_SOM") else ""
        rows.append(
            {
                "isoform": isoform,
                "source_row": source_row,
                "source_smiles": Chem.MolToSmiles(molecule),
                "connectivity": connectivity(molecule),
                "primary_som": primary_som,
                "atom_count": molecule.GetNumAtoms(),
            }
        )
external = pl.DataFrame(rows).join(challenge, on="connectivity", how="left", validate="m:1")
external.write_parquet(out / "attnsom-source-records.parquet")
overlap = external.filter(pl.col("compound_ids").is_not_null()).sort(
    "challenge_role", "compound_ids", "isoform"
)
overlap.write_parquet(out / "attnsom-challenge-overlap.parquet")
invalid_table = pl.DataFrame(invalid) if invalid else pl.DataFrame(
    schema={"isoform": pl.String, "source_row": pl.Int64}
)
invalid_table.write_parquet(out / "attnsom-invalid-records.parquet")

counts = (
    external.group_by("isoform")
    .agg(
        pl.len().alias("instances"),
        pl.col("connectivity").n_unique().alias("unique_connectivities"),
        pl.col("compound_ids").is_not_null().sum().alias("challenge_overlap_instances"),
        pl.col("overlaps_test").fill_null(False).sum().alias("test_overlap_instances"),
    )
    .sort("isoform")
)
counts.write_csv(out / "attnsom-source-counts.csv")
summary = {
    "source_commit": cfg["source_commit"],
    "instances": external.height,
    "unique_connectivities": external["connectivity"].n_unique(),
    "invalid_records": invalid_table.height,
    "challenge_overlap_instances": overlap.height,
    "challenge_overlap_connectivities": overlap["connectivity"].n_unique(),
    "challenge_train_overlap_connectivities": overlap.filter(
        pl.col("overlaps_train")
    )["connectivity"].n_unique(),
    "challenge_test_overlap_connectivities": overlap.filter(
        pl.col("overlaps_test")
    )["connectivity"].n_unique(),
    "training_policy": (
        "Exclude every challenge connectivity from ATTNSOM fitting before using its "
        "probabilities or embeddings in challenge validation."
    ),
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(counts.write_csv(), flush=True)
print(json.dumps(summary, indent=2), flush=True)
