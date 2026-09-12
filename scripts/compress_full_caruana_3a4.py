"""Replace the full-data CYP3A4 Caruana prediction with its top-K weighted models."""

import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
endpoint = cfg["endpoint"]

original = pl.read_parquet(inp / cfg["predictions"])
membership = pl.read_parquet(inp / cfg["membership"])
compounds = pl.read_parquet(inp / cfg["compounds"]).sort("compound_id")
test_mapping = pl.read_parquet(inp / cfg["test_mapping"])
test_ids = test_mapping.select("compound_id").unique().sort("compound_id")
assert test_ids.height == 750

selection = json.loads((inp / cfg["selection"]).read_text())
ranked = sorted(selection, key=lambda item: (-item["weight"], item["name"]))
retained = ranked[: cfg["k"]]
assert [item["name"] for item in retained] == cfg["expected_candidates"]
weight_sum = sum(item["weight"] for item in retained)

ids = compounds.select("compound_id")
class_table = pl.read_parquet(inp / cfg["class_emb"])
class_columns = sorted(c for c in class_table.columns if c != "compound_id")
class_features = (
    ids.join(class_table, on="compound_id", validate="1:1", maintain_order="left")
    .select(class_columns)
    .to_numpy()
    .astype(np.float32)
)
fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
morgan_features = np.stack(
    [
        fpgen.GetFingerprintAsNumPy(Chem.MolFromSmiles(smiles))
        for smiles in compounds["standardized_smiles"]
    ]
)
views = {"class": class_features, "morgan": morgan_features}
index = {compound_id: i for i, compound_id in enumerate(compounds["compound_id"])}
itest = np.array([index[c] for c in test_ids["compound_id"]])

prediction = np.zeros(test_ids.height)
manifest = []
for item in retained:
    name = item["name"]
    weight = item["weight"] / weight_sum
    booster = lgb.Booster(model_file=str(inp / cfg["models"][name]))
    p = booster.predict(views[item["view"]][itest])
    assert np.isfinite(p).all()
    prediction += weight * p
    manifest.append(
        {
            "name": name,
            "source_weight": item["weight"],
            "normalized_weight": weight,
            "view": item["view"],
        }
    )

replacement = test_ids.with_columns(
    pl.lit(endpoint).alias("endpoint"),
    pl.Series("compressed_prediction", prediction),
)
compressed = (
    original.join(replacement, on=["compound_id", "endpoint"], how="left")
    .with_columns(
        pl.coalesce("compressed_prediction", "prediction").alias("prediction")
    )
    .drop("compressed_prediction")
    .sort("endpoint", "compound_id")
)
assert compressed.height == original.height
assert compressed["prediction"].is_finite().all()
compressed.write_parquet(out / "test-predictions.parquet")
membership.write_parquet(out / "membership.parquet")
(out / "summary.json").write_text(
    json.dumps(
        {
            "endpoint": endpoint,
            "k": cfg["k"],
            "original_distinct_models": len(selection),
            "retained": manifest,
            "method": "Renormalized weighted sum of the top-K full-data Caruana models.",
            "other_endpoint_predictions": (
                "Copied unchanged from the original full-data Caruana component."
            ),
        },
        indent=2,
    )
)
