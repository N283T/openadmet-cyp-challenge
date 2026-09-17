"""Curate pinned OpenADMET metabolism datasets into the external assay corpus."""

import hashlib
import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import polars as pl
from rdkit import Chem, RDLogger

from cyp_challenge.preprocessing import standardize

EXPANSION_SPECS = {
    "LogD": ("ADME_EXPANSION_LOGD", "logD", "identity"),
    "KSOL": ("ADME_EXPANSION_KSOL", "log10_KSOL", "log10"),
    "HLM CLint": ("ADME_EXPANSION_HLM_CLINT", "log10_CLint", "log10"),
    "RLM CLint": ("ADME_EXPANSION_RLM_CLINT", "log10_CLint", "log10"),
    "MLM CLint": ("ADME_EXPANSION_MLM_CLINT", "log10_CLint", "log10"),
    "Caco-2 Permeability Papp A>B": (
        "ADME_EXPANSION_CACO2_PAPP",
        "log10_Papp_AB",
        "log10",
    ),
    "Caco-2 Permeability Efflux": (
        "ADME_EXPANSION_CACO2_EFFLUX",
        "log10_efflux_ratio",
        "log10",
    ),
    "MPPB": ("ADME_EXPANSION_MPPB", "log10_MPPB", "log10"),
    "MBPB": ("ADME_EXPANSION_MBPB", "log10_MBPB", "log10"),
    "MGMB": ("ADME_EXPANSION_MGMB", "log10_MGMB", "log10"),
}

PXR_SPECS = {
    "pEC50": ("PXR_PEC50", "pEC50", "-log10_molar"),
    "Emax_estimate (log2FC vs. baseline)": (
        "PXR_EMAX_LOG2FC",
        "Emax_log2FC",
        "identity",
    ),
    "Emax.vs.pos.ctrl_estimate (dimensionless)": (
        "PXR_EMAX_NORMALIZED",
        "Emax_normalized",
        "identity",
    ),
}


def structure_record(smiles: str) -> dict:
    RDLogger.DisableLog("rdApp.*")
    result = standardize(smiles)
    result["raw_smiles"] = smiles
    result["connectivity"] = None
    if result["compound_id"]:
        mol = Chem.MolFromSmiles(result["standardized_smiles"])
        inchikey = Chem.MolToInchiKey(mol) if mol is not None else None
        result["connectivity"] = inchikey.split("-")[0] if inchikey else None
    return result


def read_input(root: Path, name: str, *, parquet: bool = False) -> pl.DataFrame:
    exact = root / name
    paths = [exact] if exact.exists() else [
        path for path in root.iterdir() if path.name.partition("-")[2] == name
    ]
    if len(paths) != 1:
        raise ValueError(f"Expected one input named {name}; got {paths}")
    if parquet:
        return pl.read_parquet(paths[0])
    separator = "\t" if paths[0].suffix == ".tsv" else ","
    return pl.read_csv(paths[0], separator=separator, infer_schema_length=10000)


def exact_observations(
    rows: pl.DataFrame,
    *,
    source: str,
    row_key: str,
    value_column: str,
    isoform: str,
    measurement_type: str,
    assay_id: str,
    transform: str = "identity",
) -> pl.DataFrame:
    value = pl.col(value_column).cast(pl.Float64, strict=False)
    if transform == "log10":
        value = pl.when(value > 0).then(value.log10())
    return (
        rows.with_columns(value.alias("transformed_value"))
        .filter(pl.col("transformed_value").is_finite() & pl.col("connectivity").is_not_null())
        .select(
            "compound_id",
            "connectivity",
            pl.lit(isoform).alias("isoform"),
            "challenge_overlap",
            (
                pl.lit(source + ":")
                + pl.col(row_key).cast(pl.String)
                + pl.lit(":" + assay_id)
            ).alias("observation_id"),
            pl.lit(source).alias("source"),
            pl.lit(source + ":" + assay_id).alias("assay_id"),
            pl.lit(measurement_type).alias("measurement_type"),
            pl.col("transformed_value").alias("value"),
            pl.lit("=").alias("relation"),
            pl.col("transformed_value").alias("lower_bound"),
            pl.col("transformed_value").alias("upper_bound"),
            pl.lit(None, dtype=pl.Float64).alias("binary_label"),
            pl.lit(
                "Pinned primary-source OpenADMET release; challenge connectivity excluded"
            ).alias("qc"),
            pl.lit(value_column).alias("raw_column"),
            pl.col(value_column).cast(pl.String).alias("raw_value"),
            pl.lit(transform).alias("transform"),
            pl.lit(1, dtype=pl.Int32).alias("transform_direction"),
        )
    )


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

    base_compounds = read_input(inp, "external-compounds.parquet", parquet=True)
    base_observations = read_input(inp, "observations.parquet", parquet=True)
    forbidden = set(
        read_input(inp, "challenge-connectivity.parquet", parquet=True)["connectivity"].drop_nulls()
    )
    inhibition = read_input(inp, "octant-inhibition.tsv")
    reactivity = read_input(inp, "octant-reactivity.tsv")
    mass_spec = read_input(inp, "octant-mass-spec.tsv")
    expansion = read_input(inp, "expansionrx-raw.csv")
    pxr_frames = [
        read_input(inp, "pxr-train.csv"),
        read_input(inp, "pxr-phase1.csv"),
        read_input(inp, "pxr-phase2.csv"),
    ]
    pxr = pl.concat(pxr_frames, how="diagonal_relaxed").with_row_index("source_row")

    reactivity = reactivity.join(
        mass_spec.select("ocnt_batch", "standardized_smiles").unique("ocnt_batch"),
        on="ocnt_batch",
        how="left",
        validate="m:1",
    )
    source_smiles = pl.concat(
        [
            inhibition.select(pl.col("standardized_smiles").alias("raw_smiles")),
            reactivity.select(pl.col("standardized_smiles").alias("raw_smiles")),
            expansion.select(pl.col("SMILES").alias("raw_smiles")),
            pxr.select(pl.col("SMILES").alias("raw_smiles")),
        ]
    ).drop_nulls().unique()
    with ProcessPoolExecutor(max_workers=cfg["workers"]) as executor:
        records = list(executor.map(structure_record, source_smiles["raw_smiles"], chunksize=100))
    mapping = pl.DataFrame(records).unique("raw_smiles")
    mapping.write_parquet(out / "structure-provenance.parquet")

    def mapped(frame: pl.DataFrame, smiles_column: str) -> pl.DataFrame:
        return frame.join(
            mapping.select(
                pl.col("raw_smiles").alias(smiles_column),
                "compound_id",
                "standardized_smiles",
                "connectivity",
                "status",
            ),
            on=smiles_column,
            how="left",
            validate="m:1",
        ).with_columns(
            pl.col("connectivity")
            .is_in(list(forbidden))
            .fill_null(False)
            .alias("challenge_overlap")
        )

    inhibition = mapped(inhibition.rename({"standardized_smiles": "raw_smiles"}), "raw_smiles")
    reactivity = mapped(reactivity.rename({"standardized_smiles": "raw_smiles"}), "raw_smiles")
    expansion = mapped(expansion.rename({"SMILES": "raw_smiles"}), "raw_smiles")
    pxr = mapped(pxr.rename({"SMILES": "raw_smiles"}), "raw_smiles")

    observations = [
        exact_observations(
            inhibition,
            source="openadmet_octant",
            row_key="ocnt_batch",
            value_column="CYP3A4_pIC50",
            isoform="CYP3A4",
            measurement_type="pIC50",
            assay_id="CYP3A4_INHIBITION",
        )
    ]
    for enzyme in ["CYP3A4", "CYP2J2"]:
        enzyme_rows = reactivity.filter(pl.col("enzyme") == enzyme)
        for value_column, suffix in [("log2fc", "LOG2FC"), ("pct_remaining", "PCT_REMAINING")]:
            observations.append(
                exact_observations(
                    enzyme_rows,
                    source="openadmet_octant",
                    row_key="ocnt_batch",
                    value_column=value_column,
                    isoform=enzyme,
                    measurement_type=value_column,
                    assay_id=f"{enzyme}_REACTIVITY_{suffix}",
                )
            )

    # ExpansionRx supersedes/extends the older Biogen release. To avoid giving repeated
    # structures extra weight, retain only connectivities absent from the old Biogen source.
    biogen_connectivity = set(
        base_observations.filter(pl.col("source") == "biogen2023_adme")["connectivity"].drop_nulls()
    )
    expansion = expansion.with_columns(
        pl.col("connectivity").is_in(list(biogen_connectivity)).alias("older_biogen_overlap")
    )
    expansion_new = expansion.filter(~pl.col("older_biogen_overlap"))
    for raw_column, (endpoint, measurement_type, transform) in EXPANSION_SPECS.items():
        relation = pl.col(raw_column).str.extract(r"^\s*(<=|>=|<|>|=)", 1).fill_null("=")
        raw_value = pl.col(raw_column).str.replace(r"^\s*(<=|>=|<|>|=)\s*", "")
        prepared = expansion_new.with_columns(
            relation.alias("parsed_relation"),
            raw_value.cast(pl.Float64, strict=False).alias("parsed_value"),
        )
        transformed = pl.col("parsed_value")
        if transform == "log10":
            transformed = pl.when(transformed > 0).then(transformed.log10())
        prepared = prepared.with_columns(transformed.alias("transformed_value")).filter(
            pl.col("transformed_value").is_finite() & pl.col("connectivity").is_not_null()
        )
        observations.append(
            prepared.select(
                "compound_id",
                "connectivity",
                pl.lit(endpoint).alias("isoform"),
                "challenge_overlap",
                (
                    pl.lit("openadmet_expansionrx:")
                    + pl.col("Molecule Name")
                    + pl.lit(":" + endpoint)
                ).alias("observation_id"),
                pl.lit("openadmet_expansionrx").alias("source"),
                pl.lit("openadmet_expansionrx:" + endpoint).alias("assay_id"),
                pl.lit(measurement_type).alias("measurement_type"),
                pl.when(pl.col("parsed_relation") == "=")
                .then(pl.col("transformed_value"))
                .alias("value"),
                pl.col("parsed_relation").alias("relation"),
                pl.when(pl.col("parsed_relation").is_in(["=", ">", ">="]))
                .then(pl.col("transformed_value"))
                .alias("lower_bound"),
                pl.when(pl.col("parsed_relation").is_in(["=", "<", "<="]))
                .then(pl.col("transformed_value"))
                .alias("upper_bound"),
                pl.lit(None, dtype=pl.Float64).alias("binary_label"),
                pl.lit(
                    "Official ExpansionRx release; older Biogen connectivity deduplicated"
                ).alias("qc"),
                pl.lit(raw_column).alias("raw_column"),
                pl.col(raw_column).cast(pl.String).alias("raw_value"),
                pl.lit(transform).alias("transform"),
                pl.lit(1, dtype=pl.Int32).alias("transform_direction"),
            )
        )

    pxr = pxr.with_columns(
        (
            pl.col("Molecule Name")
            + pl.lit(":")
            + pl.col("source_row").cast(pl.String)
        ).alias("pxr_row_key")
    )
    for value_column, (endpoint, measurement_type, _transform) in PXR_SPECS.items():
        observations.append(
            exact_observations(
                pxr,
                source="openadmet_pxr",
                row_key="pxr_row_key",
                value_column=value_column,
                isoform=endpoint,
                measurement_type=measurement_type,
                assay_id=endpoint,
                transform="identity",
            )
        )

    all_new = pl.concat(observations, how="diagonal_relaxed")
    all_new.write_parquet(out / "all-new-observations.parquet")
    overlap = all_new.filter(pl.col("challenge_overlap"))
    accepted = all_new.filter(~pl.col("challenge_overlap"))
    overlap.write_parquet(out / "challenge-overlap-observations.parquet")
    accepted.write_parquet(out / "openadmet-observations.parquet")

    compounds = (
        pl.concat(
            [
                inhibition,
                reactivity,
                expansion_new,
                pxr,
            ],
            how="diagonal_relaxed",
        )
        .filter(pl.col("connectivity").is_not_null() & ~pl.col("challenge_overlap"))
        .select("compound_id", "standardized_smiles", "connectivity")
        .unique()
    )
    # Only retain structures that contribute at least one accepted observation.
    compounds = compounds.join(
        accepted.select("compound_id").unique(), on="compound_id", how="semi"
    )
    merged_compounds = (
        pl.concat([base_compounds, compounds]).unique("compound_id").sort("compound_id")
    )
    merged_observations = pl.concat([base_observations, accepted], how="diagonal_relaxed").sort(
        "observation_id"
    )
    assert merged_compounds["compound_id"].n_unique() == merged_compounds.height
    assert merged_observations["observation_id"].n_unique() == merged_observations.height
    assert set(merged_compounds["connectivity"]).isdisjoint(forbidden)
    assert merged_observations.join(merged_compounds, on="compound_id", how="anti").height == 0
    merged_compounds.write_parquet(out / "external-compounds.parquet")
    merged_observations.write_parquet(out / "observations.parquet")

    inventory = accepted.group_by("source", "isoform", "assay_id", "measurement_type").agg(
        pl.len().alias("observations"),
        pl.col("compound_id").n_unique().alias("compounds"),
        pl.col("value").is_not_null().sum().alias("exact"),
        pl.col("relation").is_in(["<", ">", "<=", ">="]).sum().alias("censored"),
        pl.col("value").min().alias("minimum"),
        pl.col("value").max().alias("maximum"),
    ).sort("source", "isoform")
    inventory.write_csv(out / "inventory.csv")
    source_files = {}
    for path in inp.iterdir():
        alias = path.name if path.name in cfg["source_revisions"] else path.name.partition("-")[2]
        if alias in cfg["source_revisions"]:
            source_files[alias] = hashlib.sha256(path.read_bytes()).hexdigest()
    summary = {
        "base_run": cfg["base_run"],
        "source_revisions": cfg["source_revisions"],
        "source_sha256": source_files,
        "standardized_raw_smiles": source_smiles.height,
        "standardization_failures": mapping["compound_id"].null_count(),
        "expansion_older_biogen_overlap_compounds": expansion.filter(
            pl.col("older_biogen_overlap")
        )["Molecule Name"].n_unique(),
        "challenge_overlap_observations": overlap.height,
        "new_compounds": merged_compounds.height - base_compounds.height,
        "new_observations": accepted.height,
        "merged_compounds": merged_compounds.height,
        "merged_observations": merged_observations.height,
        "heads": inventory.to_dicts(),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
