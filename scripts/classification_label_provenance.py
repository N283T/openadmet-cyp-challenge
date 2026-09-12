"""Independent raw boolean provenance, transparent rule checks and official LB comparisons."""

import hashlib
import json
import os
from pathlib import Path

import httpx
import numpy as np
import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(key):
    return pl.read_parquet(inp / cfg[key])


aliases, standard, canonical, conflicts, compounds = [
    read(k) for k in ["aliases", "standard", "labels", "conflicts", "compounds"]
]
eps = ["CYP2D6_is_TDI", "CYP3A4_is_TDI"]
canonical = canonical.filter(pl.col("endpoint").is_in(eps))
assert (
    canonical.height == 5081
    and not canonical.select("compound_id", "endpoint").is_duplicated().any()
)
assert conflicts.filter(pl.col("endpoint").is_in(eps)).is_empty()
assert all(
    hashlib.sha256(s.encode()).hexdigest() == c
    for c, s in compounds.select("compound_id", "standardized_smiles").iter_rows()
)
identity = (
    aliases.select("source_file", "molecule_name", "raw_smiles", "compound_id")
    .unique()
    .with_columns(pl.lit(True).alias("alias_present"))
)
assert (
    identity.group_by("source_file", "molecule_name", "raw_smiles")
    .agg(pl.col("compound_id").n_unique())["compound_id"]
    .max()
    == 1
)
raw_frames, stats, rule_frames = [], [], []
for filename in cfg["raw"]:
    path = inp / filename
    source = filename.split("-", 1)[1]
    table = pl.read_csv(
        path, infer_schema_length=None, null_values=["", "NaN", "nan"]
    ).with_row_index("source_row", offset=1)
    keys = table.select(
        pl.lit(source).alias("source_file"),
        pl.col("Molecule_Name").alias("molecule_name"),
        pl.col("SMILES").alias("raw_smiles"),
        "source_row",
    )
    mapped = keys.join(
        identity,
        on=["source_file", "molecule_name", "raw_smiles"],
        how="left",
        validate="m:1",
        maintain_order="left",
    )
    assert mapped.height == table.height and mapped["alias_present"].null_count() == 0
    for ep in eps:
        if ep not in table.columns:
            continue
        value = table[ep].cast(pl.String).str.to_lowercase()
        assert set(value.drop_nulls()).issubset({"true", "false", "0", "1"})
        numeric = value.replace_strict(
            {"true": 1, "false": 0, "1": 1, "0": 0}, default=None, return_dtype=pl.Int8
        )
        frame = mapped.with_columns(pl.lit(ep).alias("endpoint"), numeric.alias("raw_label"))
        raw_frames.append(frame)
        stats.append(
            dict(
                source_file=source,
                endpoint=ep,
                rows=table.height,
                labeled=int(numeric.is_not_null().sum()),
                positive=int(numeric.fill_null(0).sum()),
                missing=numeric.null_count(),
            )
        )
        iso = ep.split("_")[0]
        direct = iso + "_pIC50_direct_inhibition"
        tdi = iso + "_pIC50_TDI_condition"
        if direct in table.columns and tdi in table.columns:
            rule = frame.with_columns(
                table[direct].cast(pl.Float64).alias("direct"),
                table[tdi].cast(pl.Float64).alias("tdi"),
            )
            rule = rule.with_columns((pl.col("tdi") - pl.col("direct")).alias("shift"))
            rule = rule.with_columns(
                pl.when(pl.col("direct") > 4)
                .then(
                    pl.when(pl.col("shift") > np.log10(2))
                    .then(pl.lit("measured_positive"))
                    .otherwise(pl.lit("measured_negative"))
                )
                .when((pl.col("direct") < 4) & (pl.col("tdi") > 4 + np.log10(2)))
                .then(pl.lit("inferred_positive"))
                .when((pl.col("direct") < 4) & (pl.col("tdi") < 4))
                .then(pl.lit("assigned_negative"))
                .otherwise(pl.lit("undefined_by_simple_rule"))
                .alias("simple_category")
            )
            # Both potency values must be present; missing TDI cannot imply a negative call.
            rule = rule.with_columns(
                pl.when(
                    pl.col("direct").is_finite()
                    & pl.col("tdi").is_finite()
                    & (pl.col("simple_category") != "undefined_by_simple_rule")
                )
                .then(
                    pl.col("simple_category")
                    .is_in(["measured_positive", "inferred_positive"])
                    .cast(pl.Int8)
                )
                .otherwise(None)
                .alias("simple_label")
            )
            rule_frames.append(rule)
raw = pl.concat(raw_frames, how="vertical_relaxed")
raw.write_parquet(out / "raw-label-map.parquet")
pl.DataFrame(stats).write_csv(out / "raw-file-counts.csv")
raw.filter(pl.col("compound_id").is_null()).write_csv(out / "quarantined-label-rows.csv")
known = raw.filter(pl.col("raw_label").is_not_null() & pl.col("compound_id").is_not_null())
unique = known.select("compound_id", "endpoint", "raw_label").unique()
assert (
    unique.group_by("compound_id", "endpoint")
    .agg(pl.col("raw_label").n_unique())["raw_label"]
    .max()
    == 1
)
joined = unique.join(
    canonical.select(
        "compound_id", "endpoint", pl.col("value").cast(pl.Int8).alias("canonical_label")
    ),
    on=["compound_id", "endpoint"],
    how="full",
    coalesce=True,
    validate="1:1",
)
assert joined.height == 5081 and joined.null_count().sum_horizontal().item() == 0
assert (joined["raw_label"] == joined["canonical_label"]).all()
joined.write_parquet(out / "canonical-reconciliation.parquet")
merged = known.group_by("compound_id", "endpoint").agg(
    pl.len().alias("raw_occurrences"),
    pl.col("raw_smiles").n_unique().alias("raw_smiles_count"),
    pl.col("molecule_name").n_unique().alias("name_count"),
)
merged.write_parquet(out / "deduplication-groups.parquet")
rules = pl.concat(rule_frames, how="vertical_relaxed").unique(
    subset=["compound_id", "endpoint", "direct", "tdi", "raw_label", "simple_category"]
)
rules.write_parquet(out / "simple-rule-comparison.parquet")
rule_summary = (
    rules.group_by("endpoint", "simple_category", "raw_label", "simple_label")
    .len()
    .sort("endpoint", "simple_category")
)
rule_summary.write_csv(out / "simple-rule-summary.csv")
rules.filter(
    pl.col("raw_label").is_not_null()
    & pl.col("simple_label").is_not_null()
    & (pl.col("raw_label") != pl.col("simple_label"))
).write_csv(out / "rule-disagreements.csv")
# Source definition and public comparison snapshots remain distinct from our interpretation.
lb = []
for ep, filename in cfg["leaderboards"].items():
    tab = pl.read_csv(inp / filename, infer_schema_length=None)
    selected = tab.filter(
        pl.col("Username").is_in(["N283T", "TabICL-baseline"]) | (pl.col("Rank") == 1)
    )
    assert selected.height == 3
    lb.append(selected.with_columns(pl.lit(ep).alias("endpoint")))
pl.concat(lb, how="diagonal_relaxed").write_csv(out / "official-baseline-comparison.csv")
url = "https://raw.githubusercontent.com/OpenADMET/CYP-Challenge-Tutorial/858ae63ce79934113bccdb7fc65467de5f7b1935/evaluation/utils.py"
response = httpx.get(url, timeout=60, follow_redirects=True)
response.raise_for_status()
(out / "official-evaluation-utils.py").write_text(response.text)
summary = dict(
    raw_labeled_occurrences=known.height,
    unique_canonical_labels=joined.height,
    label_mismatches=0,
    label_conflicts=0,
    unmapped_source_rows=0,
    merged_multi_smiles_groups=merged.filter(pl.col("raw_smiles_count") > 1).height,
    merged_multi_name_groups=merged.filter(pl.col("name_count") > 1).height,
    source_files=stats,
    rule_disagreements=rules.filter(
        pl.col("raw_label").is_not_null()
        & pl.col("simple_label").is_not_null()
        & (pl.col("raw_label") != pl.col("simple_label"))
    ).height,
    utils_source_url=url,
    utils_sha256=hashlib.sha256(response.content).hexdigest(),
)
(out / "summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary), flush=True)
