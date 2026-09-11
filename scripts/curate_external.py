"""Curate source-separated CYP auxiliary readouts and exclude challenge identities."""

import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem, RDLogger

from cyp_challenge.preprocessing import standardize
from cyp_challenge.storage import write_json


def identity(smiles):
    RDLogger.DisableLog("rdApp.*")
    row = standardize(smiles)
    row["connectivity"] = None
    if row["compound_id"]:
        ik = Chem.MolToInchiKey(Chem.MolFromSmiles(row["standardized_smiles"]))
        row["connectivity"] = ik.split("-")[0] if ik else None
    return row


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

    def read(name):
        return pl.read_parquet(next(inp.glob(f"*-{name}.parquet")))

    chembl, pubchem, challenge = read("chembl-raw"), read("pubchem-raw"), read("compounds")
    # Use the original deposited structure, retaining SID/CID and all raw source columns.
    pubchem = pubchem.with_columns(
        pl.col("Panel Name").str.replace("p450-", "").str.to_uppercase().alias("isoform"),
        pl.col("Max_Response").cast(pl.Float64, strict=False).alias("raw_value"),
        pl.col("PUBCHEM_EXT_DATASOURCE_SMILES").alias("smiles"),
    )
    assert set(pubchem["isoform"]) == {"CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4", "CYP2C19"}
    assert pubchem.select("PUBCHEM_SID", "isoform").n_unique() == pubchem.height
    chembl = chembl.with_columns(
        (pl.lit("CYP") + pl.col("target_name").str.replace("Cytochrome P450 ", "")).alias(
            "isoform"
        ),
    )
    accepted = (
        pl.col("standard_type").is_in(cfg["chembl_allowed_types"])
        & ~pl.col("assay_src_id").is_in(cfg["chembl_excluded_source_ids"])
        & (pl.col("target_type") == "SINGLE PROTEIN")
        & pl.col("assay_type").is_in(cfg.get("chembl_assay_types", ["B", "F"]))
        & (pl.col("confidence_score") >= cfg["chembl_min_confidence"])
        & (pl.col("standard_relation") == "=")
        & (pl.col("standard_units") == "nM")
        & (pl.col("standard_value") > 0)
        & pl.col("pchembl_value").is_finite()
        & pl.col("data_validity_comment").is_null()
        & (pl.col("potential_duplicate").fill_null(0) == 0)
        & ~pl.col("assay_description")
        .str.contains(cfg["chembl_assay_exclusion_regex"])
        .fill_null(False)
    ).fill_null(False)
    chembl = chembl.with_columns(accepted.alias("passes_assay_policy"))
    chembl.group_by("isoform", "standard_type", "assay_src_id", "passes_assay_policy").agg(
        pl.len().alias("rows")
    ).sort("isoform", "standard_type", "assay_src_id").write_csv(out / "chembl-selection.csv")
    # Excluded records, including bounds and imported panel potency, remain inspectable.
    chembl.write_parquet(out / "chembl-selection-records.parquet")
    selected = chembl.filter(pl.col("passes_assay_policy"))
    assert np.allclose(
        9 - np.log10(selected["standard_value"]), selected["pchembl_value"], atol=0.015
    )
    all_smiles = sorted(
        set(selected["smiles"].drop_nulls())
        | set(pubchem["smiles"].drop_nulls())
        | set(challenge["standardized_smiles"])
    )
    structs = []
    with ProcessPoolExecutor(max_workers=cfg["standardization_workers"]) as executor:
        for i, row in enumerate(executor.map(identity, all_smiles, chunksize=100)):
            structs.append(row)
            if i % 2000 == 0:
                print(f"Standardized {i}/{len(all_smiles)}", flush=True)
    structures = pl.DataFrame(
        structs, infer_schema_length=None, schema_overrides={"error": pl.String}
    )
    structures.write_parquet(out / "external-standardization.parquet")
    mapping = structures.select(
        pl.col("raw_smiles").alias("smiles"),
        "compound_id",
        "standardized_smiles",
        "connectivity",
        "status",
    )
    challenge_map = challenge.join(
        mapping.rename(
            {
                "smiles": "standardized_smiles",
                "standardized_smiles": "restandardized_smiles",
                "compound_id": "recomputed_id",
            }
        ),
        on="standardized_smiles",
        validate="1:1",
    )
    assert challenge_map.height == challenge.height
    assert (challenge_map["compound_id"] == challenge_map["recomputed_id"]).all()
    assert challenge_map["connectivity"].null_count() == 0
    challenge_map.select("compound_id", "connectivity").write_parquet(
        out / "challenge-connectivity.parquet"
    )
    forbidden = challenge_map["connectivity"].unique()
    c = selected.join(mapping, on="smiles", how="left", validate="m:1").with_columns(
        pl.col("connectivity")
        .is_in(forbidden.implode())
        .fill_null(False)
        .alias("challenge_overlap")
    )
    p = pubchem.join(mapping, on="smiles", how="left", validate="m:1").with_columns(
        pl.col("connectivity")
        .is_in(forbidden.implode())
        .fill_null(False)
        .alias("challenge_overlap"),
        pl.col("raw_value").clip(*cfg["max_response_clip"]).alias("value"),
    )
    c.write_parquet(out / "chembl-curated-records.parquet")
    p.write_parquet(out / "pubchem-curated-records.parquet")
    accepted_c = c.filter(pl.col("connectivity").is_not_null() & ~pl.col("challenge_overlap"))
    accepted_p = p.filter(
        pl.col("connectivity").is_not_null()
        & ~pl.col("challenge_overlap")
        & pl.col("value").is_finite()
    )
    cols = ["compound_id", "standardized_smiles", "connectivity"]
    compounds = (
        pl.concat([accepted_c.select(cols), accepted_p.select(cols)]).unique().sort("compound_id")
    )
    assert compounds.height == compounds["compound_id"].n_unique()
    assert set(compounds["connectivity"]).isdisjoint(forbidden)
    ca = accepted_c.group_by("compound_id", "isoform", "assay_chembl_id").agg(
        pl.col("pchembl_value").median().alias("value"), pl.len().alias("source_rows")
    )
    ct = ca.group_by("compound_id", "isoform").agg(
        pl.col("value").median(), pl.col("source_rows").sum(), pl.len().alias("assays")
    )
    ct = ct.with_columns(
        (pl.lit("chembl_ic50_") + pl.col("isoform")).alias("endpoint"),
        pl.lit("chembl_ic50").alias("source"),
    )
    pt = accepted_p.group_by("compound_id", "isoform").agg(
        pl.col("value").median(), pl.len().alias("source_rows")
    )
    pt = pt.with_columns(
        (pl.lit("pubchem_max_response_") + pl.col("isoform")).alias("endpoint"),
        pl.lit("pubchem_max_response").alias("source"),
        pl.lit(1).alias("assays"),
    )
    targets = pl.concat([ct, pt], how="diagonal_relaxed").sort("compound_id", "endpoint")
    assert targets.select("compound_id", "endpoint").n_unique() == targets.height
    targets.write_parquet(out / "external-targets.parquet")
    compounds.write_parquet(out / "external-compounds.parquet")
    counts = (
        targets.group_by("source", "endpoint")
        .agg(
            pl.len().alias("n"),
            pl.col("value").mean().alias("mean"),
            pl.col("value").min().alias("min"),
            pl.col("value").max().alias("max"),
            pl.col("source_rows").sum().alias("source_rows"),
        )
        .sort("source", "endpoint")
    )
    counts.write_csv(out / "target-counts.csv")
    p.group_by("isoform", "PUBCHEM_ACTIVITY_OUTCOME").agg(
        pl.len().alias("rows"),
        pl.col("raw_value").is_finite().sum().alias("finite_response"),
        (pl.col("raw_value") != pl.col("value")).sum().alias("clipped"),
    ).sort("isoform", "PUBCHEM_ACTIVITY_OUTCOME").write_csv(out / "pubchem-response-coverage.csv")
    chembl.filter(pl.col("assay_src_id") == 7).group_by(
        "src_assay_id", "isoform", "standard_type"
    ).agg(pl.len().alias("rows")).sort("src_assay_id", "isoform").write_csv(
        out / "chembl-pubchem-assays.csv"
    )
    summary = dict(
        config=cfg,
        compounds=compounds.height,
        labels=targets.height,
        targets=counts.to_dicts(),
        chembl_accepted_rows_before_structure=c.height,
        chembl_overlap_rows=c.filter(pl.col("challenge_overlap")).height,
        pubchem_overlap_rows=p.filter(pl.col("challenge_overlap")).height,
        chembl_structure_failed_rows=c.filter(pl.col("connectivity").is_null()).height,
        pubchem_structure_failed_rows=p.filter(pl.col("connectivity").is_null()).height,
        pubchem_missing_response_rows=p.filter(
            ~pl.col("raw_value").is_finite().fill_null(False)
        ).height,
        external_sources_shared_compounds=len(set(ct["compound_id"]) & set(pt["compound_id"])),
        overlap_policy="All challenge connectivity blocks excluded, beyond fold-only exclusion",
        bounds_policy="Preserved in source records, excluded from exact-target experiment",
    )
    write_json(out / "summary.json", summary)
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
