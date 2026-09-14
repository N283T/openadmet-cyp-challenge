"""Audit exact structure overlap between prior PXR labels and CYP compounds."""

import json
import os
from pathlib import Path

import polars as pl
from rdkit import Chem, RDLogger

from cyp_challenge.preprocessing import standardize

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
RDLogger.DisableLog("rdApp.*")


def connectivity(smiles: str) -> str | None:
    molecule = Chem.MolFromSmiles(smiles)
    return Chem.MolToInchiKey(molecule).split("-")[0] if molecule else None


def source_rows(
    table: pl.DataFrame,
    source: str,
    smiles_column: str,
    value_column: str,
) -> pl.DataFrame:
    rows = []
    for row in table.select(smiles_column, value_column).iter_rows(named=True):
        normalized = standardize(row[smiles_column])
        smiles = normalized["standardized_smiles"]
        rows.append(
            {
                "source": source,
                "connectivity": connectivity(smiles) if smiles else None,
                "value": row[value_column],
                "status": normalized["status"],
            }
        )
    return pl.DataFrame(rows)


compounds = pl.read_parquet(inp / cfg["compounds"]).select("compound_id", "standardized_smiles")
test_ids = set(pl.read_parquet(inp / cfg["test_mapping"])["compound_id"])
challenge = compounds.with_columns(
    pl.col("standardized_smiles")
    .map_elements(connectivity, return_dtype=pl.String)
    .alias("connectivity"),
    pl.col("compound_id")
    .map_elements(lambda value: "test" if value in test_ids else "train", return_dtype=pl.String)
    .alias("cyp_split"),
)
sources = pl.concat(
    [
        source_rows(
            pl.read_parquet(inp / cfg["pxr_train_activity"]),
            "pxr_activity_train",
            "smiles",
            "pec50",
        ),
        source_rows(
            pl.read_parquet(inp / cfg["pxr_test_activity"]),
            "pxr_activity_released_test",
            "SMILES",
            "pEC50",
        ),
        source_rows(
            pl.read_parquet(inp / cfg["pxr_single_concentration"]),
            "pxr_single_concentration",
            "SMILES",
            "log2_fc_estimate",
        ),
    ],
    how="vertical_relaxed",
)
invalid = sources.filter(pl.col("connectivity").is_null())
valid = sources.filter(pl.col("connectivity").is_not_null())
overlap = valid.join(challenge, on="connectivity", how="inner").sort(
    "source", "cyp_split", "compound_id"
)
overlap.write_parquet(out / "overlap-records.parquet")
summary_rows = []
for source in valid["source"].unique().sort():
    source_table = valid.filter(pl.col("source") == source)
    source_overlap = overlap.filter(pl.col("source") == source)
    summary_rows.append(
        {
            "source": source,
            "raw_rows": source_table.height,
            "unique_connectivities": source_table["connectivity"].n_unique(),
            "invalid_rows": invalid.filter(pl.col("source") == source).height,
            "overlap_train_connectivities": source_overlap.filter(pl.col("cyp_split") == "train")[
                "connectivity"
            ].n_unique(),
            "overlap_test_connectivities": source_overlap.filter(pl.col("cyp_split") == "test")[
                "connectivity"
            ].n_unique(),
        }
    )
summary = pl.DataFrame(summary_rows)
summary.write_csv(out / "summary.csv")
split_conflicts = (
    challenge.group_by("connectivity")
    .agg(
        pl.col("compound_id").n_unique().alias("cyp_compounds"),
        pl.col("cyp_split").n_unique().alias("split_count"),
    )
    .filter(pl.col("split_count") > 1)
)
split_conflicts.write_parquet(out / "challenge-connectivity-split-conflicts.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "policy": (
                "Classical PXR auxiliary models will exclude every CYP connectivity. "
                "Frozen checkpoints require separate treatment if their training source overlaps."
            ),
            "challenge_connectivities_in_both_train_and_test": split_conflicts.height,
            "rows": summary_rows,
        },
        indent=2,
    )
    + "\n"
)
print(summary, flush=True)
