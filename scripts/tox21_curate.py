"""Standardize Tox21 identities and retain source-specific label alternatives."""

import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import polars as pl
from rdkit import Chem, RDLogger

from cyp_challenge.preprocessing import standardize


def identity(smiles):
    RDLogger.DisableLog("rdApp.*")
    row = standardize(smiles)
    row["connectivity"] = None
    if row["compound_id"]:
        key = Chem.MolToInchiKey(Chem.MolFromSmiles(row["standardized_smiles"]))
        row["connectivity"] = key.split("-")[0] if key else None
    return row


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    raw = pl.read_parquet(next(inp.glob("*-tox21-raw.parquet")))
    rep = pl.read_parquet(next(inp.glob("*-tox21-replicates.parquet")))
    challenge = pl.read_parquet(next(inp.glob("*-challenge-connectivity.parquet")))
    prior = pl.read_parquet(next(inp.glob("*-external-standardization.parquet")))
    smiles = sorted(set(raw["PUBCHEM_EXT_DATASOURCE_SMILES"].drop_nulls()))
    with ProcessPoolExecutor(max_workers=4) as pool:
        mapped = list(pool.map(identity, smiles, chunksize=100))
    structs = pl.DataFrame(mapped, infer_schema_length=None, schema_overrides={"error": pl.String})
    structs.write_parquet(out / "tox21-standardization.parquet")
    mapping = structs.select(
        pl.col("raw_smiles").alias("PUBCHEM_EXT_DATASOURCE_SMILES"),
        "compound_id",
        "standardized_smiles",
        "connectivity",
        "status",
    )
    t = raw.join(
        mapping, on="PUBCHEM_EXT_DATASOURCE_SMILES", how="left", validate="m:1"
    ).with_columns(
        pl.col("connectivity")
        .is_in(challenge["connectivity"].unique().implode())
        .fill_null(False)
        .alias("challenge_overlap"),
        pl.col("connectivity")
        .is_in(prior["connectivity"].drop_nulls().unique().implode())
        .fill_null(False)
        .alias("prior_external_overlap"),
    )
    t.write_parquet(out / "tox21-curated-records.parquet")
    r = rep.join(
        mapping, on="PUBCHEM_EXT_DATASOURCE_SMILES", how="left", validate="m:1"
    ).with_columns(
        (-pl.col("Fit_LogAC50").cast(pl.Float64, strict=False)).alias("pac50"),
        pl.col("Fit_R2").cast(pl.Float64, strict=False).alias("fit_r2"),
        pl.col("Fit_CurveClass").cast(pl.Float64, strict=False).alias("curve_class"),
    )
    r.write_parquet(out / "tox21-curated-replicates.parquet")
    # Overall outcomes combine replicates. Never assign Inconclusive to the inactive class.
    accepted = t.filter(pl.col("connectivity").is_not_null() & ~pl.col("challenge_overlap"))
    binary = accepted.filter(
        pl.col("PUBCHEM_ACTIVITY_OUTCOME").is_in(["Active", "Inactive"])
    ).select(
        "compound_id",
        "isoform",
        "aid",
        "PUBCHEM_SID",
        (pl.col("PUBCHEM_ACTIVITY_OUTCOME") == "Active").cast(pl.Float64).alias("value"),
    )
    binary = binary.group_by("compound_id", "isoform", "aid").agg(
        pl.col("value").mean(), pl.len().alias("n_sids")
    )
    binary = binary.with_columns(pl.lit("binary").alias("label_kind"))
    # Conservative potency: overall Active, replicate Inhibitor, negative quality curve.
    active = accepted.filter(pl.col("PUBCHEM_ACTIVITY_OUTCOME") == "Active").select(
        "aid", "PUBCHEM_SID"
    )
    potency = r.join(active, on=["aid", "PUBCHEM_SID"], validate="m:1").filter(
        (pl.col("Phenotype") == "Inhibitor")
        & pl.col("pac50").is_finite()
        & pl.col("curve_class").is_in([-1.1, -1.2, -2.1, -2.2])
    )
    potency.group_by("isoform", "Phenotype", "curve_class").len().write_csv(
        out / "potency-curve-counts.csv"
    )
    potency = potency.group_by("compound_id", "isoform", "aid").agg(
        pl.col("pac50").median().alias("value"), pl.col("PUBCHEM_SID").n_unique().alias("n_sids")
    )
    potency = potency.with_columns(pl.lit("pac50").alias("label_kind"))
    labels = pl.concat([binary, potency], how="diagonal_relaxed").with_columns(
        pl.lit("tox21_p450_glo").alias("source")
    )
    labels.write_parquet(out / "tox21-targets.parquet")
    accepted.select("compound_id", "standardized_smiles", "connectivity").unique().write_parquet(
        out / "tox21-compounds.parquet"
    )
    labels.group_by("isoform", "label_kind").agg(
        pl.len().alias("n"), pl.col("value").mean().alias("mean")
    ).sort("isoform", "label_kind").write_csv(out / "target-counts.csv")
    summary = dict(
        raw_sids=raw["PUBCHEM_SID"].n_unique(),
        raw_rows=raw.height,
        standardized_compounds=t["compound_id"].drop_nulls().n_unique(),
        challenge_overlap_compounds=t.filter(pl.col("challenge_overlap"))["compound_id"].n_unique(),
        accepted_compounds=accepted["compound_id"].n_unique(),
        new_vs_prior_compounds=accepted.filter(~pl.col("prior_external_overlap"))[
            "compound_id"
        ].n_unique(),
        labels=labels.height,
        potency_policy="Active Inhibitor; curve -1.1,-1.2,-2.1,-2.2; median",
        binary_policy="Active=1 Inactive=0; exclude Inconclusive; mean duplicates",
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
