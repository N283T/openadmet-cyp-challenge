"""Audit source-selection attrition before choosing the first auxiliary training dataset."""

import os
from pathlib import Path

import polars as pl


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    raw = pl.read_parquet(next(inp.glob("*-chembl-raw.parquet")))
    steps = [
        ("IC50", pl.col("standard_type") == "IC50"),
        ("non_pubchem", pl.col("assay_src_id") != 7),
        ("single_protein", pl.col("target_type") == "SINGLE PROTEIN"),
        ("B_or_F", pl.col("assay_type").is_in(["B", "F"])),
        ("confidence8", pl.col("confidence_score") >= 8),
        (
            "exact_positive_nM",
            (pl.col("standard_relation") == "=")
            & (pl.col("standard_units") == "nM")
            & (pl.col("standard_value") > 0),
        ),
        ("pchembl", pl.col("pchembl_value").is_finite()),
        ("no_validity_flag", pl.col("data_validity_comment").is_null()),
        ("not_duplicate", pl.col("potential_duplicate").fill_null(0) == 0),
    ]
    rows, f = [], raw
    for step, expression in steps:
        f = f.filter(expression)
        rows.extend(
            f.group_by("target_name")
            .agg(pl.len().alias("rows"))
            .with_columns(pl.lit(step).alias("step"))
            .to_dicts()
        )
    pl.DataFrame(rows).write_csv(out / "filter-funnel.csv")
    f.select("assay_chembl_id", "assay_description", "target_name").unique().write_csv(
        out / "eligible-assays.csv"
    )
    qhts = pl.read_parquet(next(inp.glob("*-pubchem-raw.parquet")))
    qhts.group_by("Compound QC").agg(pl.len().alias("rows")).write_csv(out / "pubchem-qc.csv")
    print(pl.DataFrame(rows), flush=True)


if __name__ == "__main__":
    main()
