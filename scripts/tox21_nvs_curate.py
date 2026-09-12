"""Link official DSSTox structures to NVS assays and prepare separate targets."""

import json
import os
from pathlib import Path

import polars as pl
from rdkit import Chem, RDLogger

from cyp_challenge.preprocessing import standardize

RDLogger.DisableLog("rdApp.*")
inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
t = pl.read_parquet(next(inp.glob("*-nvs-raw.parquet")))
forbidden = pl.read_parquet(next(inp.glob("*-challenge-connectivity.parquet")))["connectivity"]
rows = []
for mol in Chem.SDMolSupplier(str(next(inp.glob("*-INVITRODBV3_20181017_v2000.sdf")))):
    if mol is None:
        continue
    props = {k: mol.GetProp(k) for k in mol.GetPropNames()}
    cas = props.get("Substance_CASRN")
    if cas not in set(t["casn"].drop_nulls()):
        continue
    smiles = props.get("Structure_SMILES") or Chem.MolToSmiles(mol)
    row = standardize(smiles)
    key = (
        Chem.MolToInchiKey(Chem.MolFromSmiles(row["standardized_smiles"]))
        if row["compound_id"]
        else None
    )
    rows.append({**props, **row, "casn": cas, "connectivity": key.split("-")[0] if key else None})
s = pl.DataFrame(rows, infer_schema_length=None, schema_overrides={"error": pl.String})
s.write_parquet(out / "nvs-dsstox-structures.parquet")
mapping = (
    s.filter(pl.col("compound_id").is_not_null())
    .select("casn", "compound_id", "standardized_smiles", "connectivity")
    .unique()
)
ambiguous = mapping.group_by("casn").len().filter(pl.col("len") > 1)["casn"]
mapping = mapping.filter(~pl.col("casn").is_in(ambiguous.implode()))
r = t.join(mapping, on="casn", how="left", validate="m:1").with_columns(
    pl.col("connectivity")
    .is_in(forbidden.unique().implode())
    .fill_null(False)
    .alias("challenge_overlap"),
    pl.col("hitc").cast(pl.Float64, strict=False).alias("hitcall"),
    (6 - pl.col("modl_ga").cast(pl.Float64, strict=False)).alias("pac50"),
    pl.col("aenm").str.replace("NVS_ADME_h", "").alias("isoform"),
)
r.write_parquet(out / "nvs-curated-records.parquet")
accepted = r.filter(
    pl.col("compound_id").is_not_null()
    & ~pl.col("challenge_overlap")
    & (pl.col("tested_conc_unit") == "uM")
    & (pl.col("gsid_rep") == "1")
    & pl.col("hitcall").is_in([0, 1])
)
accepted.select("compound_id", "standardized_smiles", "connectivity").unique().write_parquet(
    out / "nvs-compounds.parquet"
)
keys = ["compound_id", "isoform", "aeid"]
b = (
    accepted.group_by(keys)
    .agg(pl.col("hitcall").mean().alias("value"), pl.len().alias("n_source_rows"))
    .with_columns(pl.lit("binary").alias("label_kind"))
)
p = (
    accepted.filter(
        (pl.col("hitcall") == 1)
        & pl.col("pac50").is_finite()
        & (pl.col("flag_length").cast(pl.Float64, strict=False).fill_null(0) == 0)
    )
    .group_by(keys)
    .agg(pl.col("pac50").median().alias("value"), pl.len().alias("n_source_rows"))
    .with_columns(pl.lit("pac50").alias("label_kind"))
)
labels = pl.concat([b, p]).with_columns(pl.lit("nvs_adme_v33").alias("source"))
labels.write_parquet(out / "nvs-targets.parquet")
counts = (
    labels.group_by("isoform", "label_kind")
    .agg(pl.len().alias("n"), pl.col("value").mean().alias("mean"))
    .sort("isoform", "label_kind")
)
counts.write_csv(out / "target-counts.csv")
summary = dict(
    raw_rows=t.height,
    raw_cas=t["casn"].drop_nulls().n_unique(),
    structure_matched_cas=r.filter(pl.col("compound_id").is_not_null())["casn"].n_unique(),
    ambiguous_cas=len(ambiguous),
    accepted_compounds=accepted["compound_id"].n_unique(),
    challenge_overlap_compounds=r.filter(pl.col("challenge_overlap"))["compound_id"].n_unique(),
    label_rows=labels.height,
    unit_policy="uM only; pAC50=6-modl_ga",
    selection="gsid_rep=1 and hitc in0,1; potency additionallyhitc1 finite unflagged",
)
(out / "summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary), flush=True)
print(counts, flush=True)
