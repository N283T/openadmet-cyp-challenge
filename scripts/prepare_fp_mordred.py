"""Assemble exact baseline Morgan+Mordred raw features without fitting preprocessing."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
compounds = pl.read_parquet(
    next(p for p in inp.iterdir() if p.name.partition("-")[2] == "compounds.parquet")
).sort("compound_id")
shards = sorted(inp.glob("*-mordred-*.parquet"))
assert shards
mordred = pl.concat(pl.read_parquet(path) for path in shards)
columns = [name for name in mordred.columns if name != "compound_id"]
assert len(columns) == 1613
assert mordred["compound_id"].n_unique() == mordred.height
joined = compounds.select("compound_id", "standardized_smiles").join(
    mordred, on="compound_id", how="left", validate="1:1", maintain_order="left"
)
assert set(compounds["compound_id"]) == set(mordred["compound_id"])
generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fp = np.stack(
    [
        generator.GetFingerprintAsNumPy(Chem.MolFromSmiles(smiles))
        for smiles in joined["standardized_smiles"]
    ]
)
assert fp.shape == (6896, 2048)
fp_names = [f"fp_{i:04d}" for i in range(2048)]
result = pl.DataFrame(fp, schema=fp_names).insert_column(0, joined["compound_id"])
# Keep source descriptor null/nonfinite cells; preprocessing remains fold-local.
result = result.hstack(joined.select(columns).rename({name: "mordred_" + name for name in columns}))
assert result.height == 6896 and result.width == 3662
result.write_parquet(out / "fp-mordred-features.parquet")
metadata = dict(
    config=cfg,
    rows=result.height,
    features=3661,
    fingerprint_dimensions=2048,
    descriptor_dimensions=1613,
    descriptor_source_shards=len(shards),
    fingerprint_columns=fp_names,
    descriptor_columns=["mordred_" + name for name in columns],
    original_descriptor_columns=columns,
    preprocessing=(
        "None fitted. Baseline: nonfinite->NaN; training-only finite/nonconstant filter; "
        "training median imputation of descriptors; fingerprints first."
    ),
    feature_order="Morgan bit0..2047 followed by original Mordred sourcecolumn order",
    missing_descriptor_cells=int(joined.select(columns).null_count().to_numpy().sum()),
)
(out / "summary.json").write_text(json.dumps(metadata, indent=2))
print(
    json.dumps({key: value for key, value in metadata.items() if not key.endswith("columns")}),
    flush=True,
)
