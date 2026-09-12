"""Curate six primary-author Biogen ADME measurements and merge onto the AZ corpus."""

import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import polars as pl
from adme_az_curate import identity

SPECS = [
    dict(
        raw_column="LOG HLM_CLint (mL/min/kg)",
        isoform="ADME_BIOGEN_HLM_CLINT",
        measurement_type="LOG_CL_mL_min_kg",
        source_assay_id="HLM_CLINT",
        assay_organism="Homo sapiens",
        unit="mL/min/kg",
        assay_system="human liver microsomes",
    ),
    dict(
        raw_column="LOG RLM_CLint (mL/min/kg)",
        isoform="ADME_BIOGEN_RLM_CLINT",
        measurement_type="LOG_CL_mL_min_kg",
        source_assay_id="RLM_CLINT",
        assay_organism="Rattus norvegicus",
        unit="mL/min/kg",
        assay_system="rat liver microsomes",
    ),
    dict(
        raw_column="LOG MDR1-MDCK ER (B-A/A-B)",
        isoform="ADME_BIOGEN_MDCK_EFFLUX",
        measurement_type="LOG_efflux_BA_AB_ratio",
        source_assay_id="MDR1_MDCK_ER",
        assay_organism=None,
        unit="B-A/A-B ratio",
        assay_system="MDR1-MDCK efflux",
    ),
    dict(
        raw_column="LOG SOLUBILITY PH 6.8 (ug/mL)",
        isoform="ADME_BIOGEN_SOLUBILITY68",
        measurement_type="LOG_solubility_ug_mL",
        source_assay_id="SOLUBILITY_PH68",
        assay_organism=None,
        unit="ug/mL",
        assay_system="solubility pH6.8",
    ),
    dict(
        raw_column="LOG PLASMA PROTEIN BINDING (HUMAN) (% unbound)",
        isoform="ADME_BIOGEN_PPB_HUMAN",
        measurement_type="LOG_unbound_percent",
        source_assay_id="HUMAN_PPB",
        assay_organism="Homo sapiens",
        unit="% unbound",
        assay_system="human plasma protein binding",
    ),
    dict(
        raw_column="LOG PLASMA PROTEIN BINDING (RAT) (% unbound)",
        isoform="ADME_BIOGEN_PPB_RAT",
        measurement_type="LOG_unbound_percent",
        source_assay_id="RAT_PPB",
        assay_organism="Rattus norvegicus",
        unit="% unbound",
        assay_system="rat plasma protein binding",
    ),
]


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

    def path(name):
        return next(p for p in inp.iterdir() if p.name.partition("-")[2] == name)

    def read(name):
        return pl.read_parquet(path(name + ".parquet"))

    source = json.loads(path("source.json").read_text())
    sha = next(r["sha256"] for r in source["files"] if r["file"] == "ADME_public_set_3521.csv")
    raw = pl.read_csv(path("ADME_public_set_3521.csv"), infer_schema=False)
    assert raw["Internal ID"].n_unique() == raw.height
    raw.write_parquet(out / "biogen-raw.parquet")
    mapping = read("az-structure-provenance")
    missing = sorted(set(raw["SMILES"]) - set(mapping["raw_smiles"]))
    with ProcessPoolExecutor(max_workers=cfg["workers"]) as ex:
        extra = list(ex.map(identity, missing, chunksize=100))
    if extra:
        mapping = pl.concat(
            [mapping, pl.DataFrame(extra, schema=mapping.schema)], how="vertical_relaxed"
        )
    mapping.write_parquet(out / "biogen-structure-provenance.parquet")
    r = raw.join(
        mapping.select(
            pl.col("raw_smiles").alias("SMILES"),
            "compound_id",
            "standardized_smiles",
            "connectivity",
            "status",
        ),
        on="SMILES",
        how="left",
        validate="m:1",
    )
    forbidden = read("challenge-connectivity")["connectivity"].unique()
    r = r.with_columns(
        pl.col("connectivity")
        .is_in(forbidden.implode())
        .fill_null(False)
        .alias("challenge_overlap")
    )
    r.write_parquet(out / "biogen-curated-compounds.parquet")
    specs = pl.DataFrame(SPECS)
    long = r.unpivot(
        index=[c for c in r.columns if c not in specs["raw_column"]],
        on=specs["raw_column"].to_list(),
        variable_name="raw_column",
        value_name="raw_value",
    ).join(specs, on="raw_column", validate="m:1")
    long = long.with_columns(
        pl.col("raw_value").str.extract(r"^\s*(<=|>=|<|>|=)", 1).fill_null("=").alias("relation"),
        pl.col("raw_value")
        .str.replace(r"^\s*(<=|>=|<|>|=)\s*", "")
        .cast(pl.Float64, strict=False)
        .alias("logged_value"),
    )
    long.write_parquet(out / "biogen-curated-observations.parquet")
    usable = long.filter(pl.col("connectivity").is_not_null() & pl.col("logged_value").is_finite())
    observations = usable.select(
        "compound_id",
        "connectivity",
        "isoform",
        "challenge_overlap",
        (pl.lit("biogen2023:") + pl.col("Internal ID") + ":" + pl.col("source_assay_id")).alias(
            "observation_id"
        ),
        pl.lit(cfg["source"]).alias("source"),
        (pl.lit("biogen2023:") + pl.col("source_assay_id")).alias("assay_id"),
        "measurement_type",
        "relation",
        pl.when(pl.col("relation") == "=").then(pl.col("logged_value")).alias("value"),
        pl.when(pl.col("relation").is_in(["=", ">", ">="]))
        .then(pl.col("logged_value"))
        .alias("lower_bound"),
        pl.when(pl.col("relation").is_in(["=", "<", "<="]))
        .then(pl.col("logged_value"))
        .alias("upper_bound"),
        pl.lit(None, dtype=pl.Float64).alias("binary_label"),
        pl.lit(None, dtype=pl.String).alias("hitcall"),
        pl.lit(None, dtype=pl.Float64).alias("max_response"),
        pl.lit(None, dtype=pl.Float64).alias("fit_r2"),
        pl.lit(None, dtype=pl.Float64).alias("curve_class"),
        pl.lit("Primary author numerical logged value; no inferred censoring or threshold").alias(
            "qc"
        ),
        pl.lit(False).alias("preincubation_mentioned"),
        "source_assay_id",
        "assay_organism",
        "assay_system",
        pl.col("unit").alias("original_property_unit"),
        "raw_column",
        "raw_value",
        pl.lit("identity_already_logged").alias("transform"),
        pl.lit(1).alias("transform_direction"),
        pl.lit(source["commit"]).alias("source_commit"),
        pl.lit(sha).alias("source_csv_sha256"),
    )
    observations.filter("challenge_overlap").write_parquet(
        out / "biogen-challenge-overlap-observations.parquet"
    )
    accepted = observations.filter(~pl.col("challenge_overlap"))
    accepted.write_parquet(out / "biogen-observations.parquet")
    c = (
        usable.filter(~pl.col("challenge_overlap"))
        .select("compound_id", "standardized_smiles", "connectivity")
        .unique()
    )
    c.write_parquet(out / "biogen-compounds.parquet")
    prior_c, prior_o = read("external-compounds"), read("observations")
    merged_c = pl.concat([prior_c, c]).unique().sort("compound_id")
    merged_o = pl.concat([prior_o, accepted], how="diagonal_relaxed").sort("observation_id")
    assert merged_c["compound_id"].n_unique() == merged_c.height
    assert merged_o["observation_id"].n_unique() == merged_o.height
    assert set(merged_c["connectivity"]).isdisjoint(forbidden)
    assert merged_o.join(merged_c, on="compound_id", how="anti").height == 0
    merged_c.write_parquet(out / "external-compounds.parquet")
    merged_o.write_parquet(out / "observations.parquet")
    counts = accepted.group_by(
        "isoform", "assay_id", "assay_organism", "original_property_unit"
    ).agg(
        pl.len().alias("observations"),
        pl.col("compound_id").n_unique().alias("compounds"),
        pl.col("value").is_not_null().sum().alias("continuous"),
        pl.col("relation").is_in(["<", ">", "<=", ">="]).sum().alias("censored"),
        pl.col("value").min().alias("minimum"),
        pl.col("value").max().alias("maximum"),
    )
    counts.write_csv(out / "biogen-inventory.csv")
    specs.write_json(out / "biogen-head-metadata.json")
    summary = dict(
        source=source,
        base_run=cfg["base_run"],
        raw_compounds=raw.height,
        raw_nonmissing_values=long["raw_value"].is_not_null().sum(),
        new_structure_mappings=len(missing),
        biogen_compounds=c.height,
        biogen_observations=accepted.height,
        new_compounds=merged_c.height - prior_c.height,
        overlap_observations=observations.filter("challenge_overlap").height,
        merged_compounds=merged_c.height,
        merged_observations=merged_o.height,
        merged_continuous=merged_o["value"].is_not_null().sum(),
        merged_binary=merged_o["binary_label"].is_not_null().sum(),
        censor_tokens=long.filter(pl.col("relation") != "=").height,
        continuous_rank_heads=counts.filter(pl.col("continuous") >= 5).height,
        heads=counts.to_dicts(),
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
