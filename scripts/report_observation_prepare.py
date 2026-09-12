"""Build observed direct/TDI potency, native signed Emax and TDI-screen targets."""

import json
import os
from pathlib import Path

import polars as pl
from scipy.stats import spearmanr

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def path(i):
    return next(inp.glob(f"{i:03d}-*"))


compounds, aliases, labels = [pl.read_parquet(path(i)) for i in range(3)]
mapping = aliases.select(
    pl.col("molecule_name").alias("Molecule_Name"),
    pl.col("raw_smiles").alias("SMILES"),
    "compound_id",
).unique()
rows = []
for ep in ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]:
    rows.append(
        labels.filter(pl.col("endpoint") == ep + "_pIC50_direct_inhibition").select(
            "compound_id",
            pl.lit(ep).alias("enzyme"),
            pl.lit("direct").alias("kind"),
            "value",
            "conf_low",
            "conf_high",
            pl.lit(None, dtype=pl.Float64).alias("concentration_M"),
        )
    )
for index, suffix, kind in [
    (3, "pIC50_TDI_condition", "tdi"),
    (4, "EmaxVsPosCtrl_TDI_condition", "emax"),
]:
    source = pl.read_csv(path(index), infer_schema=False).join(
        mapping, on=["Molecule_Name", "SMILES"], how="left", validate="m:1"
    )
    unmapped = source.filter(pl.col("compound_id").is_null())
    unmapped.write_csv(out / f"unmapped-{kind}.csv")
    assert unmapped.is_empty(), (kind, unmapped.select("Molecule_Name", "SMILES").head(10))
    for ep in ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]:
        head = ep + "_" + suffix
        t = source.select(
            "compound_id",
            pl.lit(ep).alias("enzyme"),
            pl.lit(kind).alias("kind"),
            pl.col(head).cast(pl.Float64, strict=False).alias("value"),
            pl.col(head + "_conf_low").cast(pl.Float64, strict=False).alias("conf_low"),
            pl.col(head + "_conf_high").cast(pl.Float64, strict=False).alias("conf_high"),
        )
        rows.append(
            t.filter(pl.col("value").is_finite())
            .group_by("compound_id", "enzyme", "kind")
            .agg(pl.col("value", "conf_low", "conf_high").median())
            .with_columns(pl.lit(None, dtype=pl.Float64).alias("concentration_M"))
        )
screen = pl.read_csv(path(5)).join(
    mapping, on=["Molecule_Name", "SMILES"], how="left", validate="m:1"
)
assert screen["compound_id"].null_count() == 0
assert (screen["concentration_M"] > 0).all()
assert (
    screen.group_by("compound_id", "enzyme")
    .agg(pl.col("concentration_M").n_unique())["concentration_M"]
    .max()
    == 1
)
rows.append(
    screen.group_by("compound_id", "enzyme")
    .agg(pl.col("log2fc_estimate").median().alias("value"), pl.col("concentration_M").first())
    .with_columns(
        pl.lit("logfc").alias("kind"),
        pl.lit(None, dtype=pl.Float64).alias("conf_low"),
        pl.lit(None, dtype=pl.Float64).alias("conf_high"),
    )
)
targets = pl.concat(rows, how="diagonal_relaxed").sort("compound_id", "enzyme", "kind")
assert not targets.select("compound_id", "enzyme", "kind").is_duplicated().any()
assert set(targets["compound_id"]).issubset(set(compounds["compound_id"]))
targets.write_parquet(out / "observation-targets.parquet")
stats = (
    targets.group_by("enzyme", "kind")
    .agg(
        pl.len().alias("n"),
        pl.col("value").min().alias("min"),
        pl.col("value").median().alias("median"),
        pl.col("value").max().alias("max"),
    )
    .sort("enzyme", "kind")
)
stats.write_csv(out / "distributions.csv")
paired = targets.pivot(on="kind", index=["compound_id", "enzyme"], values="value").join(
    screen.group_by("compound_id", "enzyme").agg(pl.col("concentration_M").first()),
    on=["compound_id", "enzyme"],
    how="left",
    validate="1:1",
)
checks = []
for (ep,), h in paired.group_by("enzyme"):
    h = h.filter(pl.all_horizontal(pl.col("tdi", "emax", "logfc").is_finite()))
    p, e, c, y = [h[k].to_numpy() for k in ["tdi", "emax", "concentration_M", "logfc"]]
    response = e / (1 + 10 ** (-p) / c)
    checks.append(dict(enzyme=ep, n=h.height, spearman=float(spearmanr(response, y).statistic)))
(out / "summary.json").write_text(
    json.dumps(
        dict(
            distributions=stats.to_dicts(),
            observed_curve_logfc_correlations=checks,
            meaning=(
                "Screen is active TDI condition; signed Emax retained. "
                "Hill response mapped by learned per-CYP positive affine scale, "
                "not treated as exact raw signal."
            ),
        ),
        indent=2,
    )
)
print(stats)
print(checks)
