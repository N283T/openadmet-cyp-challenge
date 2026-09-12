"""Compare fixed auxiliary-source replacements in the submitted OOF recipe."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
reference = pl.read_parquet(inp / cfg["reference"]).sort("compound_id", "endpoint", "fold")
keys = ["compound_id", "endpoint", "fold", "y_true", "conf_low", "conf_high", "is_neighbor"]


def values(filename):
    h = pl.read_parquet(inp / filename).sort("compound_id", "endpoint", "fold")
    assert h.height == 6384
    assert h.select(keys).equals(reference.select(keys))
    x = h["prediction"].to_numpy()
    assert np.isfinite(x).all()
    return x


old = [values(f) for f in cfg["original_components"]]
assert np.allclose(np.mean(old, axis=0), reference["prediction"].to_numpy(), atol=1e-12, rtol=0)
replacements = {int(i): values(f) for i, f in cfg["replacements"].items()}
variants = {"original": np.mean(old, axis=0)}
recipes = {"original": dict(original_indices=list(range(5)), replacement_indices=[])}
for i, new in replacements.items():
    x = old.copy()
    x[i] = new
    name = f"replace_component_{i}"
    variants[name] = np.mean(x, axis=0)
    recipes[name] = dict(original_indices=[j for j in range(5) if j != i], replacement_indices=[i])
updated = [replacements.get(i, x) for i, x in enumerate(old)]
variants["replace_all_auxiliary_components"] = np.mean(updated, axis=0)
recipes["replace_all_auxiliary_components"] = dict(replacement_indices=sorted(replacements))
if cfg.get("addition"):
    extra = values(cfg["addition"])
    variants["original_plus_primary_lgbm"] = np.mean(old + [extra], axis=0)
    variants["replaced_plus_primary_lgbm"] = np.mean(updated + [extra], axis=0)
    recipes["original_plus_primary_lgbm"] = dict(equal_components=6, updated=False)
    recipes["replaced_plus_primary_lgbm"] = dict(equal_components=6, updated=True)
spec = importlib.util.spec_from_file_location("official", inp / cfg["scorer"])
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
rows = []
for name, p in variants.items():
    h = reference.with_columns(pl.Series("prediction", p), pl.lit(name).alias("method"))
    h.write_parquet(out / f"{name}.parquet")
    for subset in ["all", "neighbors"]:
        g = h if subset == "all" else h.filter(pl.col("is_neighbor"))
        for (ep,), z in g.group_by("endpoint"):
            score = float(
                scorer.rae_soft_threshold_absolute_error(
                    z["y_true"].to_numpy(),
                    z["prediction"].to_numpy(),
                    y_true_lower=z["conf_low"].to_numpy(),
                    y_true_upper=z["conf_high"].to_numpy(),
                )
            )
            rows.append(dict(method=name, subset=subset, endpoint=ep, n=z.height, ST_RAE=score))
metrics = pl.DataFrame(rows)
metrics.write_csv(out / "metrics.csv")
macro = metrics.group_by("method", "subset").agg(pl.col("ST_RAE").mean()).sort("subset", "ST_RAE")
macro.write_csv(out / "macro.csv")
(out / "manifest.json").write_text(
    json.dumps(
        dict(
            config=cfg,
            recipes=recipes,
            exact_label_alignment=True,
            original_reconstruction=True,
            weight_search=False,
            status="development comparison; no test inference or submission",
        ),
        indent=2,
    )
)
print(macro)
