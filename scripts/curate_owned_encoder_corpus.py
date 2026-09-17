"""Build leak-excluded endpoint variants for project-owned encoder pretraining."""

import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import polars as pl
from rdkit import Chem, RDLogger

from cyp_challenge.preprocessing import standardize

CLASSIFICATION = {
    "hia_hou",
    "pgp_broccatelli",
    "bioavailability_ma",
    "cyp2c19_veith",
    "cyp2d6_veith",
    "cyp3a4_veith",
    "cyp1a2_veith",
    "cyp2c9_veith",
    "cyp2c9_substrate_carbonmangels",
    "cyp2d6_substrate_carbonmangels",
    "cyp3a4_substrate_carbonmangels",
    "bbb_martins",
    "b3db_classification",
}

CYP_TDC = {
    "cyp2c19_veith",
    "cyp2d6_veith",
    "cyp3a4_veith",
    "cyp1a2_veith",
    "cyp2c9_veith",
    "cyp2c9_substrate_carbonmangels",
    "cyp2d6_substrate_carbonmangels",
    "cyp3a4_substrate_carbonmangels",
}

# These are already present from their primary ChEMBL/AZ acquisition and are
# retained in the raw TDC snapshot only, not duplicated as extra supervision.
TDC_DUPLICATES = {
    "lipophilicity_astrazeneca",
    "clearance_hepatocyte_az",
    "clearance_microsome_az",
}


def input_path(root: Path, name: str) -> Path:
    exact = root / name
    matches = (
        [exact]
        if exact.exists()
        else [path for path in root.iterdir() if path.name.partition("-")[2] == name]
    )
    if len(matches) != 1:
        raise ValueError(f"Expected one input {name}; got {matches}")
    return matches[0]


def structure_record(smiles: str) -> dict:
    RDLogger.DisableLog("rdApp.*")
    record = standardize(smiles)
    record["raw_smiles"] = smiles
    record["connectivity"] = None
    if record["compound_id"]:
        molecule = Chem.MolFromSmiles(record["standardized_smiles"])
        key = Chem.MolToInchiKey(molecule) if molecule is not None else None
        record["connectivity"] = key.split("-")[0] if key else None
    return record


def source_isoform(name: str) -> str:
    if name.startswith("cyp"):
        token = name.split("_", 1)[0].upper().replace("CYP1A2", "CYP1A2")
        return token
    return f"ADME_TDC_{name.upper()}"


def tdc_rows(path: Path, name: str) -> pl.DataFrame:
    frame = pl.read_csv(path, separator="\t", infer_schema_length=0).with_row_index("source_row")
    smiles = "Drug" if "Drug" in frame.columns else "X"
    frame = frame.select(
        "source_row",
        pl.col(smiles).cast(pl.String).alias("raw_smiles"),
        pl.col("Y").cast(pl.Float64, strict=False).alias("raw_y"),
        *([pl.col("Species").cast(pl.String)] if "Species" in frame.columns else []),
    ).filter(pl.col("raw_smiles").is_not_null() & pl.col("raw_y").is_finite())
    assay = pl.lit(f"tdc:{name}")
    if "Species" in frame.columns:
        assay = assay + pl.lit(":") + pl.col("Species").fill_null("unspecified")
    isoform = source_isoform(name)
    is_classification = name in CLASSIFICATION
    return frame.select(
        "source_row",
        "raw_smiles",
        pl.lit(name).alias("dataset"),
        pl.lit("tdc_adme").alias("source"),
        assay.alias("assay_id"),
        pl.lit(isoform).alias("isoform"),
        pl.lit("binary" if is_classification else "continuous").alias("measurement_type"),
        pl.when(pl.lit(not is_classification)).then(pl.col("raw_y")).alias("value"),
        pl.when(pl.lit(not is_classification)).then(pl.col("raw_y")).alias("lower_bound"),
        pl.when(pl.lit(not is_classification)).then(pl.col("raw_y")).alias("upper_bound"),
        pl.when(pl.lit(is_classification)).then(pl.col("raw_y")).alias("binary_label"),
        pl.col("raw_y").cast(pl.String).alias("raw_value"),
    )


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

    base_compounds = pl.read_parquet(input_path(inp, "external-compounds.parquet"))
    base_observations = pl.read_parquet(input_path(inp, "observations.parquet"))
    forbidden = set(
        pl.read_parquet(input_path(inp, "challenge-connectivity.parquet"))[
            "connectivity"
        ].drop_nulls()
    )

    tdc = []
    for name in cfg["tdc_datasets"]:
        tdc.append(tdc_rows(input_path(inp, f"{name}.tab"), name))
    tdc = pl.concat(tdc, how="diagonal_relaxed")

    nov = pl.read_parquet(input_path(inp, "novartis-train-test.parquet")).with_row_index(
        "source_row"
    )
    nov_long = []
    for column, isoform, measurement in [
        ("log_kobs", "CYP3A4", "log_kobs"),
        ("pIC50_CYP3A4", "CYP3A4", "pIC50"),
        ("pIC50_CYP2C9", "CYP2C9", "pIC50"),
        ("pIC50_CYP2D6", "CYP2D6", "pIC50"),
    ]:
        nov_long.append(
            nov.select(
                "source_row",
                pl.col("Structure").cast(pl.String).alias("raw_smiles"),
                pl.lit("novartis_cyp3a4").alias("dataset"),
                pl.lit("polaris_novartis_cyp").alias("source"),
                pl.lit(f"polaris_novartis_cyp:{column}").alias("assay_id"),
                pl.lit(isoform).alias("isoform"),
                pl.lit(measurement).alias("measurement_type"),
                pl.col(column).cast(pl.Float64, strict=False).alias("value"),
                pl.col(column).cast(pl.Float64, strict=False).alias("lower_bound"),
                pl.col(column).cast(pl.Float64, strict=False).alias("upper_bound"),
                pl.lit(None, dtype=pl.Float64).alias("binary_label"),
                pl.col(column).cast(pl.String).alias("raw_value"),
            ).filter(pl.col("raw_smiles").is_not_null() & pl.col("value").is_finite())
        )
    new_rows = pl.concat([tdc, *nov_long], how="diagonal_relaxed")

    raw_structures = new_rows.select("raw_smiles").unique()
    with ProcessPoolExecutor(max_workers=cfg["workers"]) as pool:
        records = list(pool.map(structure_record, raw_structures["raw_smiles"], chunksize=100))
    structures = pl.DataFrame(records).unique("raw_smiles")
    structures.write_parquet(out / "structure-provenance.parquet")
    new_rows = new_rows.join(
        structures.select(
            "raw_smiles", "compound_id", "standardized_smiles", "connectivity", "status"
        ),
        on="raw_smiles",
        how="left",
        validate="m:1",
    ).with_columns(
        pl.col("connectivity").is_in(list(forbidden)).fill_null(False).alias("challenge_overlap")
    )
    accepted = new_rows.filter(
        pl.col("compound_id").is_not_null() & ~pl.col("challenge_overlap")
    ).with_columns(
        (
            pl.col("source")
            + pl.lit(":")
            + pl.col("dataset")
            + pl.lit(":")
            + pl.col("source_row").cast(pl.String)
            + pl.lit(":")
            + pl.col("measurement_type")
        ).alias("observation_id"),
        pl.lit("=").alias("relation"),
        pl.lit("Pinned public endpoint; challenge connectivity excluded").alias("qc"),
    )
    accepted.select(
        "compound_id",
        "connectivity",
        "isoform",
        "challenge_overlap",
        "observation_id",
        "source",
        "assay_id",
        "measurement_type",
        "value",
        "relation",
        "lower_bound",
        "upper_bound",
        "binary_label",
        "qc",
        "raw_value",
    ).write_parquet(out / "new-observations.parquet")
    new_rows.filter(pl.col("challenge_overlap")).write_parquet(
        out / "challenge-overlap-observations.parquet"
    )

    new_compounds = accepted.select("compound_id", "standardized_smiles", "connectivity").unique(
        "compound_id"
    )
    compounds = (
        pl.concat([base_compounds, new_compounds], how="vertical_relaxed")
        .unique("compound_id", keep="first")
        .sort("compound_id")
    )
    compounds.write_parquet(out / "external-compounds.parquet")

    new_obs = pl.read_parquet(out / "new-observations.parquet")
    all_assay = pl.concat([base_observations, new_obs], how="diagonal_relaxed")
    all_assay = all_assay.filter(
        ~(
            (pl.col("source") == "tdc_adme")
            & pl.col("assay_id").str.extract(r"^tdc:([^:]+)", 1).is_in(sorted(TDC_DUPLICATES))
        )
    )
    all_assay.write_parquet(out / "observations-assay.parquet")

    compact = all_assay.with_columns(
        pl.when(pl.col("source") == "chembl36")
        .then(
            pl.lit("chembl36_compact:")
            + pl.col("isoform")
            + pl.lit(":")
            + pl.col("measurement_type")
        )
        .otherwise(pl.col("assay_id"))
        .alias("assay_id")
    )
    compact.write_parquet(out / "observations-compact-all.parquet")

    cyp_sources = ["chembl36", "pubchem1851", "tox21_p450_glo", "tdc_adme", "polaris_novartis_cyp"]
    compact_cyp = compact.filter(
        pl.col("source").is_in(cyp_sources)
        & (
            (pl.col("source") != "tdc_adme")
            | pl.col("assay_id").str.extract(r"^tdc:([^:]+)", 1).is_in(sorted(CYP_TDC))
        )
    )
    compact_cyp.write_parquet(out / "observations-compact-cyp.parquet")

    endpoint_only = compact.filter(pl.col("source") != "chembl36")
    endpoint_only.write_parquet(out / "observations-endpoint-only.parquet")

    inventories = []
    for name, frame in [
        ("assay", all_assay),
        ("compact_all", compact),
        ("compact_cyp", compact_cyp),
        ("endpoint_only", endpoint_only),
    ]:
        inventories.append(
            frame.group_by("source")
            .agg(
                pl.len().alias("rows"),
                pl.col("compound_id").n_unique().alias("compounds"),
                pl.col("assay_id").n_unique().alias("assays"),
            )
            .with_columns(pl.lit(name).alias("variant"))
        )
    pl.concat(inventories).sort("variant", "source").write_csv(out / "inventory.csv")
    summary = {
        "compounds": compounds.height,
        "new_accepted_rows": new_obs.height,
        "new_overlap_rows": new_rows.filter(pl.col("challenge_overlap")).height,
        "variants": {
            "assay": all_assay.height,
            "compact_all": compact.height,
            "compact_cyp": compact_cyp.height,
            "endpoint_only": endpoint_only.height,
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
