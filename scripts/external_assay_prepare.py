"""Build assay-aware external observations with continuous, censored and binary labels."""

import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import polars as pl
from rdkit import Chem, RDLogger

from cyp_challenge.preprocessing import standardize


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

    c, p = read("chembl-raw"), read("pubchem-curated-records")
    mapping = read(cfg.get("structure_mapping", "external-standardization"))
    forbidden = read("challenge-connectivity")["connectivity"].unique()
    eligible = (
        (pl.col("target_type") == "SINGLE PROTEIN")
        & (pl.col("target_organism") == "Homo sapiens")
        & (pl.col("assay_organism").is_null() | (pl.col("assay_organism") == "Homo sapiens"))
        & (pl.col("confidence_score") >= cfg["chembl_min_confidence"])
        & pl.col("assay_type").is_in(["A", "B", "F"])
        & pl.col("standard_type").is_in(cfg["chembl_types"])
        & ~pl.col("assay_src_id").is_in(cfg["chembl_exclude_sources"])
        & pl.col("data_validity_comment").is_null()
        & (pl.col("potential_duplicate").fill_null(0) == 0)
        & ~pl.col("assay_description").str.contains(cfg["tdi_regex"]).fill_null(False)
    ).fill_null(False)
    c = c.with_columns(eligible.alias("passes_assay_policy"))
    c.write_parquet(out / "chembl-selection-records.parquet")
    selected = c.filter("passes_assay_policy")
    missing = sorted(set(selected["smiles"].drop_nulls()) - set(mapping["raw_smiles"]))
    extra = []
    with ProcessPoolExecutor(max_workers=cfg["standardization_workers"]) as ex:
        for i, r in enumerate(ex.map(identity, missing, chunksize=100)):
            extra.append(r)
            if i % 2000 == 0:
                print(f"Standardized {i}/{len(missing)} new structures", flush=True)
    if extra:
        mapping = pl.concat(
            [mapping, pl.DataFrame(extra, schema=mapping.schema)], how="vertical_relaxed"
        )
    mapping.write_parquet(out / "structure-provenance.parquet")
    c = selected.join(
        mapping.select(
            pl.col("raw_smiles").alias("smiles"),
            "compound_id",
            "standardized_smiles",
            "connectivity",
            "status",
        ),
        on="smiles",
        validate="m:1",
        how="left",
    )
    c = c.with_columns(
        (pl.lit("CYP") + pl.col("target_name").str.replace("Cytochrome P450 ", "")).alias(
            "isoform"
        ),
        pl.col("connectivity")
        .is_in(forbidden.implode())
        .fill_null(False)
        .alias("challenge_overlap"),
        pl.col("assay_description")
        .str.contains("(?i)pre[- ]?incubat")
        .fill_null(False)
        .alias("preincubation_mentioned"),
    )
    c.write_parquet(out / "chembl-curated-records.parquet")
    p.write_parquet(out / "pubchem-curated-records.parquet")
    potency = (
        (pl.col("standard_units") == "nM")
        & (pl.col("standard_value") > 0)
        & pl.col("standard_value").is_finite()
    )
    c = c.with_columns(
        pl.when(potency).then(9 - pl.col("standard_value").log10()).alias("potency_p"),
        pl.col("standard_relation")
        .replace({"<": ">", ">": "<", "<=": ">=", ">=": "<="})
        .alias("relation"),
    )
    c = c.with_columns(
        pl.when(pl.col("relation") == "=").then(pl.col("potency_p")).alias("value"),
        pl.when(pl.col("relation").is_in(["=", ">", ">="]))
        .then(pl.col("potency_p"))
        .alias("lower_bound"),
        pl.when(pl.col("relation").is_in(["=", "<", "<="]))
        .then(pl.col("potency_p"))
        .alias("upper_bound"),
    )
    threshold = cfg["classification_threshold_p"]
    # Do not collapse unknown censored potencies to their bound. Binary only if unambiguous.
    c = c.with_columns(
        pl.when(pl.col("lower_bound") >= threshold)
        .then(1.0)
        .when(pl.col("upper_bound") < threshold)
        .then(0.0)
        .when((pl.col("relation") == "<") & (pl.col("upper_bound") == threshold))
        .then(0.0)
        .when((pl.col("relation") == "=") & (pl.col("value") == threshold))
        .then(1.0)
        .when(
            pl.col("potency_p").is_null()
            & pl.col("activity_comment").str.to_lowercase().is_in(["inactive", "not active"])
        )
        .then(0.0)
        .otherwise(None)
        .alias("binary_label"),
    )
    common = ["compound_id", "isoform", "connectivity", "challenge_overlap"]
    co = c.select(
        *common,
        (pl.lit("chembl:") + pl.col("activity_id").cast(pl.String)).alias("observation_id"),
        pl.lit("chembl36").alias("source"),
        (
            pl.lit("chembl36:")
            + pl.col("assay_chembl_id")
            + ":"
            + pl.col("isoform")
            + ":"
            + pl.col("standard_type")
        ).alias("assay_id"),
        pl.col("standard_type").alias("measurement_type"),
        "value",
        "relation",
        "lower_bound",
        "upper_bound",
        "binary_label",
        pl.col("activity_comment").alias("hitcall"),
        pl.lit(None, dtype=pl.Float64).alias("max_response"),
        pl.lit(None, dtype=pl.Float64).alias("fit_r2"),
        pl.lit(None, dtype=pl.Float64).alias("curve_class"),
        pl.col("data_validity_comment").alias("qc"),
        "preincubation_mentioned",
    )
    p = p.with_columns(
        pl.when(pl.col("PUBCHEM_ACTIVITY_OUTCOME") != "Inactive")
        .then(-pl.col("Fit_LogAC50").cast(pl.Float64, strict=False))
        .alias("pac50")
    )
    po = p.select(
        *common,
        (pl.lit("pubchem:1851:") + pl.col("PUBCHEM_SID") + ":" + pl.col("Panel ID")).alias(
            "observation_id"
        ),
        pl.lit("pubchem1851").alias("source"),
        (pl.lit("pubchem1851:") + pl.col("Panel ID") + ":" + pl.col("isoform") + ":AC50").alias(
            "assay_id"
        ),
        pl.lit("AC50").alias("measurement_type"),
        pl.when(pl.col("pac50").is_finite()).then(pl.col("pac50")).alias("value"),
        pl.when(pl.col("pac50").is_finite()).then(pl.lit("=")).alias("relation"),
        pl.when(pl.col("pac50").is_finite()).then(pl.col("pac50")).alias("lower_bound"),
        pl.when(pl.col("pac50").is_finite()).then(pl.col("pac50")).alias("upper_bound"),
        pl.when(pl.col("PUBCHEM_ACTIVITY_OUTCOME") == "Active")
        .then(1.0)
        .when(pl.col("PUBCHEM_ACTIVITY_OUTCOME") == "Inactive")
        .then(0.0)
        .alias("binary_label"),
        pl.col("PUBCHEM_ACTIVITY_OUTCOME").alias("hitcall"),
        pl.col("Max_Response").cast(pl.Float64, strict=False).alias("max_response"),
        pl.col("Fit_R2").cast(pl.Float64, strict=False).alias("fit_r2"),
        pl.col("Fit_CurveClass").cast(pl.Float64, strict=False).alias("curve_class"),
        pl.col("Compound QC").alias("qc"),
        pl.lit(False).alias("preincubation_mentioned"),
    )
    obs = pl.concat([co, po], how="vertical_relaxed").sort("observation_id")
    assert obs["observation_id"].n_unique() == obs.height
    obs.write_parquet(out / "all-observations.parquet")
    usable = obs.filter(
        pl.col("connectivity").is_not_null()
        & (
            pl.col("value").is_not_null()
            | pl.col("binary_label").is_not_null()
            | pl.col("lower_bound").is_not_null()
            | pl.col("upper_bound").is_not_null()
            | pl.col("max_response").is_not_null()
        )
    )
    overlaps = usable.filter("challenge_overlap")
    accepted = usable.filter(~pl.col("challenge_overlap"))
    accepted.write_parquet(out / "observations.parquet")
    overlaps.write_parquet(out / "challenge-overlap-observations.parquet")
    structs = pl.concat(
        [
            c.select("compound_id", "standardized_smiles", "connectivity"),
            p.select("compound_id", "standardized_smiles", "connectivity"),
        ]
    ).unique()
    compounds = structs.join(
        accepted.select("compound_id").unique(), on="compound_id", how="semi"
    ).sort("compound_id")
    assert compounds["compound_id"].n_unique() == compounds.height
    assert set(compounds["connectivity"]).isdisjoint(forbidden)
    compounds.write_parquet(out / "external-compounds.parquet")
    inventory = (
        accepted.group_by("source", "isoform", "measurement_type")
        .agg(
            pl.len().alias("observations"),
            pl.col("compound_id").n_unique().alias("compounds"),
            pl.col("assay_id").n_unique().alias("assays"),
            pl.col("value").is_not_null().sum().alias("continuous"),
            pl.col("binary_label").is_not_null().sum().alias("binary"),
            (pl.col("binary_label") == 0).sum().alias("inactive"),
            (pl.col("binary_label") == 1).sum().alias("active"),
            (pl.col("relation").is_in(["<", ">", "<=", ">="])).sum().alias("censored"),
        )
        .sort("source", "isoform", "measurement_type")
    )
    inventory.write_csv(out / "inventory.csv")
    assays = accepted.group_by("assay_id", "source", "isoform", "measurement_type").agg(
        pl.len().alias("observations"),
        pl.col("compound_id").n_unique().alias("compounds"),
        pl.col("value").is_not_null().sum().alias("continuous"),
        pl.col("value").drop_nulls().n_unique().alias("distinct_values"),
        pl.col("binary_label").is_not_null().sum().alias("binary"),
    )
    assays.write_parquet(out / "assay-inventory.parquet")
    summary = dict(
        config=cfg,
        compounds=compounds.height,
        observations=accepted.height,
        assays=assays.height,
        overlap_observations=overlaps.height,
        new_standardized_smiles=len(missing),
        ranking_assays_at_least_two_values=assays.filter(
            (pl.col("continuous") >= 2) & (pl.col("distinct_values") >= 2)
        ).height,
        inventory=inventory.to_dicts(),
        schema={k: str(v) for k, v in accepted.schema.items()},
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
