"""Create source ablations while retaining the identical molecular split universe."""

import hashlib
import json
import os
import shutil
from pathlib import Path

import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
compound_path = next(inp.glob("*-external-compounds.parquet"))
c = pl.read_parquet(compound_path)
o = pl.read_parquet(next(inp.glob("*-observations.parquet")))
assert c.height == c["compound_id"].n_unique() == 35320
assert set(cfg["sources"]) <= set(o["source"])
selected = o.filter(pl.col("source").is_in(cfg["sources"]))
assert selected.height == selected["observation_id"].n_unique()
assert not selected["challenge_overlap"].any()
assert selected.join(c.select("compound_id"), on="compound_id", how="anti").height == 0
shutil.copyfile(compound_path, out / "external-compounds.parquet")
source_digest = hashlib.sha256(compound_path.read_bytes()).hexdigest()
assert (
    hashlib.sha256((out / "external-compounds.parquet").read_bytes()).hexdigest() == source_digest
)
selected.write_parquet(out / "observations.parquet")
usable = selected.filter(pl.col("value").is_not_null() | pl.col("binary_label").is_not_null())
counts = selected.group_by("source", "isoform", "measurement_type").agg(
    pl.len().alias("observations"),
    pl.col("compound_id").n_unique().alias("observed_compounds"),
    pl.col("assay_id").n_unique().alias("assays"),
    pl.col("value").is_not_null().sum().alias("continuous"),
    pl.col("binary_label").is_not_null().sum().alias("binary"),
)
counts.sort("source", "isoform", "measurement_type").write_csv(out / "inventory.csv")
summary = dict(
    config=cfg,
    compounds=c.height,
    compound_pool_sha256=source_digest,
    source_compound_pool_identical=True,
    observations=selected.height,
    observed_compounds=selected["compound_id"].n_unique(),
    labelled_compounds=usable["compound_id"].n_unique(),
    unused_pool_compounds=c.height - usable["compound_id"].n_unique(),
    continuous=selected["value"].is_not_null().sum(),
    binary=selected["binary_label"].is_not_null().sum(),
    binary_compounds=selected.filter(pl.col("binary_label").is_not_null())[
        "compound_id"
    ].n_unique(),
    assays=selected["assay_id"].n_unique(),
    label_policy="Preserve masks, relations and quality flags; filter sources only",
    split_policy="Full35320compoundpool byte-identical to external3; same split universe",
    inventory=counts.to_dicts(),
)
(out / "summary.json").write_text(json.dumps(summary, indent=2))
(out / "README.md").write_text(
    "# External source ablation\n\n"
    f"Sources: {', '.join(cfg['sources'])}.\n\n"
    f"The full {c.height:,}-compound pool is byte-identical to the external3 input. "
    f"Only {summary['labelled_compounds']:,} compounds have usable labels in this view; "
    f"{summary['unused_pool_compounds']:,} pool compounds are unused by these objectives. "
    "Keep the same trainer seed and connectivity split.\n"
)
print(json.dumps({k: v for k, v in summary.items() if k != "inventory"}), flush=True)
