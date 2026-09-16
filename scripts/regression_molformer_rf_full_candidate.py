"""Full-fit the selected MoLFormer RF and replace only D6 in the live candidate."""

import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor

from cyp_challenge.tracks import validate

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def path(index: int) -> Path:
    matches = list(inp.glob(f"{index:03d}-*"))
    if len(matches) != 1:
        raise ValueError(f"Expected one input at {index}, found {matches}")
    return matches[0]


measurements = pl.read_parquet(path(0))
scores = pl.read_parquet(path(1))
features = pl.read_parquet(path(2))
auxiliary = pl.read_parquet(path(3))
aliases = pl.read_parquet(path(4))
current = pl.read_csv(path(5)).sort("Molecule_Name")
reference = path(6)

endpoint = cfg["endpoint"]
score_columns = sorted(column for column in scores.columns if column.startswith("score_"))
feature_columns = sorted(
    column for column in features.columns if column.startswith("molformer_cls_mean_")
)
aux_columns = [f"aux_{index}" for index in range(4)]
columns = score_columns + feature_columns + aux_columns
if (len(score_columns), len(feature_columns), len(columns)) != (47, 1536, 1587):
    raise ValueError("Unexpected deployment feature dimensions")

static = scores.join(features, on="compound_id", validate="1:1")
train_aux = auxiliary.filter(pl.col("role") == "cross_fitted_train")
test_aux = auxiliary.filter(pl.col("role") == "full_fit_test")
train = (
    measurements.filter(pl.col("endpoint") == endpoint)
    .select("compound_id", "value")
    .join(static, on="compound_id", validate="1:1")
    .join(train_aux.select("compound_id", *aux_columns), on="compound_id", validate="1:1")
    .sort("compound_id")
)
test_ids = (
    aliases.filter(pl.col("role") == "test")
    .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
    .unique()
)
test = (
    test_ids.join(static, on="compound_id", validate="1:1")
    .join(test_aux.select("compound_id", *aux_columns), on="compound_id", validate="1:1")
    .sort("Molecule_Name")
)
if train.height != cfg["expected_train_rows"] or test.height != 750:
    raise ValueError(f"Unexpected train/test rows: {train.height}/{test.height}")
x_train = train.select(columns).to_numpy().astype(np.float32)
x_test = test.select(columns).to_numpy().astype(np.float32)
y_train = train["value"].to_numpy()
if not np.isfinite(x_train).all() or not np.isfinite(x_test).all():
    raise ValueError("Non-finite deployment features")

model = RandomForestRegressor(
    n_estimators=cfg["tree_estimators"],
    min_samples_leaf=cfg["min_samples_leaf"],
    max_features=cfg["max_features"],
    n_jobs=cfg["threads"],
    random_state=cfg["seed"],
)
model.fit(x_train, y_train)
raw = model.predict(x_test)
joblib.dump({"model": model, "columns": columns}, out / "model.joblib")

old = current[endpoint].to_numpy()
matched = (raw - raw.mean()) / raw.std() * old.std() + old.mean()
blended = (1 - cfg["blend_weight"]) * old + cfg["blend_weight"] * matched
candidate = current.with_columns(pl.Series(endpoint, blended))
candidate.write_parquet(out / "predictions.parquet")
candidate.write_csv(out / "predictions.csv")
validation = validate(out / "predictions.parquet", "regression", reference)
if not validation["valid"]:
    raise ValueError(validation)

rank_old = pl.Series(old).rank("average").to_numpy()
rank_new = pl.Series(blended).rank("average").to_numpy()
test.select("compound_id", "Molecule_Name").with_columns(
    pl.Series("rf_prediction", raw),
    pl.Series("matched_rf_prediction", matched),
    pl.Series("old_d6_prediction", old),
    pl.Series("new_d6_prediction", blended),
    pl.Series("old_rank", rank_old),
    pl.Series("new_rank", rank_new),
).write_parquet(out / "d6-test-predictions.parquet")

summary = {
    "config": cfg,
    "validation": validation,
    "train_rows": train.height,
    "test_rows": test.height,
    "feature_dimensions": len(columns),
    "old_mean": float(old.mean()),
    "old_std": float(old.std()),
    "raw_mean": float(raw.mean()),
    "raw_std": float(raw.std()),
    "new_mean": float(blended.mean()),
    "new_std": float(blended.std()),
    "old_new_spearman": float(spearmanr(old, blended).statistic),
    "old_rf_spearman": float(spearmanr(old, raw).statistic),
    "changed_ranks": int(np.sum(rank_old != rank_new)),
    "mean_absolute_change": float(np.abs(blended - old).mean()),
    "maximum_absolute_change": float(np.abs(blended - old).max()),
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary), flush=True)
