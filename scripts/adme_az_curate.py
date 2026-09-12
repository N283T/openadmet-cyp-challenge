"""Curate measured ADME properties as separate heads and merge with frozen CYP corpus."""

import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import polars as pl
from rdkit import Chem, RDLogger

from cyp_challenge.preprocessing import standardize

SPECS = {
    "CHEMBL3301363": dict(
        endpoint="ADME_LOGD74",
        measurement_type="LogD7.4",
        units=None,
        species=None,
        transform="identity",
        direction=1,
    ),
    "CHEMBL3301364": dict(
        endpoint="ADME_SOLUBILITY74",
        measurement_type="log10_solubility_M",
        units="nM",
        species=None,
        transform="log10_nM_minus9",
        direction=1,
    ),
    "CHEMBL3301365": dict(
        endpoint="ADME_PPB_HUMAN",
        measurement_type="log10_unbound_percent",
        units="%",
        species="Homo sapiens",
        transform="log10_100_minus_bound_percent",
        direction=-1,
    ),
    "CHEMBL3301370": dict(
        endpoint="ADME_HLM_CLINT",
        measurement_type="log10_CL_mL_min_g",
        units="mL.min-1.g-1",
        species="Homo sapiens",
        transform="log10",
        direction=1,
    ),
    "CHEMBL3301372": dict(
        endpoint="ADME_HEPATOCYTE_CLINT",
        measurement_type="log10_CL_uL_min_million_cells",
        units="uL.min-1.(10^6cells)-1",
        species="Homo sapiens",
        transform="log10",
        direction=1,
    ),
}


def identity(smiles):
    RDLogger.DisableLog("rdApp.*")
    r = standardize(smiles)
    r["connectivity"] = None
    if r["compound_id"]:
        ik = Chem.MolToInchiKey(Chem.MolFromSmiles(r["standardized_smiles"]))
        r["connectivity"] = ik.split("-")[0] if ik else None
    return r


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

    def read(name):
        return pl.read_parquet(
            next(p for p in inp.iterdir() if p.name.partition("-")[2] == name + ".parquet")
        )

    raw = read("az-raw")
    mapping = read("structure-provenance")
    missing = sorted(set(raw["raw_smiles"].drop_nulls()) - set(mapping["raw_smiles"]))
    with ProcessPoolExecutor(max_workers=cfg["workers"]) as ex:
        extra = list(ex.map(identity, missing, chunksize=100))
    if extra:
        mapping = pl.concat(
            [mapping, pl.DataFrame(extra, schema=mapping.schema)], how="vertical_relaxed"
        )
    mapping.write_parquet(out / "az-structure-provenance.parquet")
    r = raw.join(
        mapping.select(
            "raw_smiles", "compound_id", "standardized_smiles", "connectivity", "status"
        ),
        on="raw_smiles",
        validate="m:1",
        how="left",
    )
    forbidden = read("challenge-connectivity")["connectivity"].unique()
    r = r.with_columns(
        pl.col("connectivity")
        .is_in(forbidden.implode())
        .fill_null(False)
        .alias("challenge_overlap")
    )
    rows = []
    for assay, spec in SPECS.items():
        g = r.filter(pl.col("source_assay_id") == assay)
        assert g.height and g["standard_units"].unique().to_list() == [spec["units"]]
        assert g["assay_organism"].unique().to_list() == [spec["species"]]
        g = g.with_columns(
            (
                pl.col("standard_value").is_finite()
                & pl.col("data_validity_comment").is_null()
                & (pl.col("potential_duplicate").fill_null(0) == 0)
            ).alias("quality_accepted")
        )
        v = pl.col("standard_value")
        if spec["transform"] == "identity":
            transformed = v
        elif spec["transform"] == "log10_nM_minus9":
            transformed = pl.when(v > 0).then(v.log10() - 9)
        elif spec["transform"] == "log10_100_minus_bound_percent":
            transformed = pl.when((v >= 0) & (v < 100)).then((100 - v).log10())
        else:
            transformed = pl.when(v > 0).then(v.log10())
        relation = pl.col("standard_relation").fill_null("=")
        if spec["direction"] == -1:
            relation = relation.replace({"<": ">", ">": "<", "<=": ">=", ">=": "<="})
        g = g.with_columns(
            transformed.alias("transformed_value"),
            relation.alias("relation"),
            pl.lit(spec["endpoint"]).alias("isoform"),
            pl.lit(spec["measurement_type"]).alias("measurement_type"),
            pl.lit(spec["transform"]).alias("transform"),
            pl.lit(spec["direction"]).alias("transform_direction"),
            pl.col("standard_relation").is_null().alias("null_relation_assumed_exact"),
        )
        rows.append(g)
    records = pl.concat(rows, how="vertical_relaxed")
    records.write_parquet(out / "az-curated-records.parquet")
    usable = records.filter(
        pl.col("quality_accepted")
        & pl.col("connectivity").is_not_null()
        & pl.col("transformed_value").is_finite()
        & pl.col("relation").is_in(["=", "<", ">", "<=", ">="])
    )
    observations = usable.select(
        "compound_id",
        "connectivity",
        "isoform",
        "challenge_overlap",
        (pl.lit("chembl36:") + pl.col("activity_id").cast(pl.String)).alias("observation_id"),
        pl.lit(cfg["source_name"]).alias("source"),
        (pl.lit("chembl36_az:") + pl.col("source_assay_id") + ":" + pl.col("isoform")).alias(
            "assay_id"
        ),
        "measurement_type",
        "relation",
        pl.when(pl.col("relation") == "=").then(pl.col("transformed_value")).alias("value"),
        pl.when(pl.col("relation").is_in(["=", ">", ">="]))
        .then(pl.col("transformed_value"))
        .alias("lower_bound"),
        pl.when(pl.col("relation").is_in(["=", "<", "<="]))
        .then(pl.col("transformed_value"))
        .alias("upper_bound"),
        pl.lit(None, dtype=pl.Float64).alias("binary_label"),
        pl.lit(None, dtype=pl.String).alias("hitcall"),
        pl.lit(None, dtype=pl.Float64).alias("max_response"),
        pl.lit(None, dtype=pl.Float64).alias("fit_r2"),
        pl.lit(None, dtype=pl.Float64).alias("curve_class"),
        pl.lit("source27 allowlisted ADME assay; unflagged; original metadata retained").alias(
            "qc"
        ),
        pl.lit(False).alias("preincubation_mentioned"),
        "source_assay_id",
        "standard_units",
        "standard_relation",
        "standard_value",
        "assay_organism",
        "transform",
        "transform_direction",
        "null_relation_assumed_exact",
    )
    observations.filter("challenge_overlap").write_parquet(
        out / "az-challenge-overlap-observations.parquet"
    )
    accepted = observations.filter(~pl.col("challenge_overlap"))
    accepted.write_parquet(out / "az-observations.parquet")
    c = (
        usable.filter(~pl.col("challenge_overlap"))
        .select("compound_id", "standardized_smiles", "connectivity")
        .unique()
    )
    c.write_parquet(out / "az-compounds.parquet")
    prior_c, prior_o = read("external-compounds"), read("observations")
    merged_c = pl.concat([prior_c, c]).unique().sort("compound_id")
    merged_o = pl.concat([prior_o, accepted], how="diagonal_relaxed").sort("observation_id")
    assert merged_c["compound_id"].n_unique() == merged_c.height
    assert merged_o["observation_id"].n_unique() == merged_o.height
    assert set(merged_c["connectivity"]).isdisjoint(forbidden)
    assert merged_o.join(merged_c, on="compound_id", how="anti").height == 0
    assert accepted["binary_label"].null_count() == accepted.height
    merged_c.write_parquet(out / "external-compounds.parquet")
    merged_o.write_parquet(out / "observations.parquet")
    inventory = accepted.group_by(
        "isoform", "measurement_type", "source_assay_id", "standard_units", "assay_organism"
    ).agg(
        pl.len().alias("observations"),
        pl.col("compound_id").n_unique().alias("compounds"),
        pl.col("value").is_not_null().sum().alias("continuous"),
        pl.col("relation").is_in(["<", ">", "<=", ">="]).sum().alias("censored"),
        pl.col("value").min().alias("minimum"),
        pl.col("value").max().alias("maximum"),
    )
    inventory.write_csv(out / "az-inventory.csv")
    (out / "head-metadata.json").write_text(json.dumps(SPECS, indent=2))
    summary = dict(
        raw_rows=raw.height,
        raw_molecules=raw["molregno"].n_unique(),
        standardized_smiles_added=len(missing),
        az_compounds=c.height,
        az_observations=accepted.height,
        az_continuous=accepted["value"].is_not_null().sum(),
        az_censored=accepted["relation"].is_in(["<", ">", "<=", ">="]).sum(),
        challenge_overlap_rows=observations.filter("challenge_overlap").height,
        new_compounds=merged_c.height - prior_c.height,
        merged_compounds=merged_c.height,
        merged_observations=merged_o.height,
        merged_continuous=merged_o["value"].is_not_null().sum(),
        merged_binary=merged_o["binary_label"].is_not_null().sum(),
        base_run=cfg["merge_base_run"],
        heads=inventory.to_dicts(),
        usable_continuous_ranking_heads=inventory.filter(pl.col("continuous") >= 5).height,
        binary_tasks_added=0,
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
