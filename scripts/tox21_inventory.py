"""Inventory raw Tox21 assays and produce replicate-preserving tables."""

import json
import os
import re
from pathlib import Path

import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
desc = json.loads(next(inp.glob("*-assay-descriptions.json")).read_text())
print("description keys", list(desc), flush=True)
metadata = []
for item in desc["PC_AssayContainer"]:
    d = item["assay"]["descr"]
    aid = d["aid"]["id"]
    print(aid, d["name"], flush=True)
    metadata.append(
        {
            "aid": aid,
            "name": d["name"],
            "description": d.get("description"),
            "protocol": d.get("protocol"),
            "results": d.get("results"),
        }
    )
(out / "assay-metadata.json").write_text(json.dumps(metadata, indent=2))
allrows, reps, counts = [], [], []
for m in metadata:
    aid = m["aid"]
    isoform = (
        re.search(r"CYP\s*(1A2|2C19|2C9|2D6|3A4)", m["name"], re.I).group().upper().replace(" ", "")
    )
    t = pl.read_csv(next(inp.glob(f"*-aid{aid}.csv")), infer_schema=False).filter(
        pl.col("PUBCHEM_SID").cast(pl.Int64, strict=False).is_not_null()
    )
    t = t.with_columns(pl.lit(aid).alias("aid"), pl.lit(isoform).alias("isoform"))
    allrows.append(t)
    print(aid, isoform, t.height, "unique SID", t["PUBCHEM_SID"].n_unique(), flush=True)
    for col in ["PUBCHEM_ACTIVITY_OUTCOME", "Compound QC-Replicate_1", "Phenotype-Replicate_1"]:
        if col in t.columns:
            print(t.group_by(col).len().to_dicts(), flush=True)
    for rep in [1, 2, 3]:
        suffix = f"-Replicate_{rep}"
        cols = [c for c in t.columns if c.endswith(suffix)]
        if not cols:
            continue
        r = t.select(
            "aid",
            "isoform",
            "PUBCHEM_SID",
            "PUBCHEM_CID",
            "PUBCHEM_EXT_DATASOURCE_SMILES",
            "PUBCHEM_ACTIVITY_OUTCOME",
            *cols,
        ).rename({c: c.removesuffix(suffix) for c in cols})
        r = r.with_columns(pl.lit(rep).alias("replicate"))
        reps.append(r)
    counts.extend(t.group_by("aid", "isoform", "PUBCHEM_ACTIVITY_OUTCOME").len().to_dicts())
pl.concat(allrows, how="diagonal_relaxed").write_parquet(out / "tox21-raw.parquet")
pl.concat(reps, how="diagonal_relaxed").write_parquet(out / "tox21-replicates.parquet")
pl.DataFrame(counts).write_csv(out / "counts.csv")
