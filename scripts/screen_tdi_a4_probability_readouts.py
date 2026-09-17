"""Cross-fit CYP3A4 pIC50 readouts from leakage-safe TDI probability features."""

import importlib.util
import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from catboost import CatBoostRegressor
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

spec = importlib.util.spec_from_file_location("official", inp / "scorer.py")
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)

base = (
    pl.read_parquet(inp / "regression-oof.parquet")
    .filter(pl.col("endpoint") == cfg["endpoint"])
    .sort("fold", "compound_id")
)
tdi_selected = (
    pl.read_parquet(inp / "tdi-oof.parquet")
    .filter(pl.col("model") == cfg["tdi_model"])
    .select("compound_id", pl.col("prediction").alias("tdi_equal7"))
)
members = pl.read_parquet(inp / "tdi-member-oof.parquet").drop("fold", "y_true")
member_columns = [column for column in members.columns if column != "compound_id"]
data = base.join(tdi_selected, on="compound_id", validate="1:1").join(
    members, on="compound_id", validate="1:1"
)
assert data.height == base.height == 2300
assert np.isfinite(data.select("tdi_equal7", *member_columns).to_numpy()).all()

methods = ["scalar_linear", "scalar_isotonic"]
methods += [f"member_ridge_{alpha:g}" for alpha in cfg["ridge_alphas"]]
methods += ["member_lgbm", "member_catboost", "member_extratrees"]
rows = {method: [] for method in methods}
fits = []
y = data["y_true"].to_numpy()
folds = data["fold"].to_numpy()
x_scalar = data.select("tdi_equal7").to_numpy()
x_members = data.select(member_columns).to_numpy()

for fold in range(5):
    valid = folds == fold
    train = ~valid
    models = {
        "scalar_linear": Ridge(alpha=1.0),
        "scalar_isotonic": IsotonicRegression(out_of_bounds="clip"),
        **{
            f"member_ridge_{alpha:g}": make_pipeline(StandardScaler(), Ridge(alpha=alpha))
            for alpha in cfg["ridge_alphas"]
        },
        "member_lgbm": lgb.LGBMRegressor(
            objective="regression_l1",
            n_estimators=300,
            num_leaves=7,
            min_child_samples=30,
            learning_rate=0.03,
            reg_lambda=3.0,
            n_jobs=cfg["num_threads"],
            random_state=42 + fold,
            verbosity=-1,
            deterministic=True,
            force_col_wise=True,
        ),
        "member_catboost": CatBoostRegressor(
            iterations=400,
            depth=4,
            learning_rate=0.03,
            l2_leaf_reg=5,
            loss_function="MAE",
            random_seed=42 + fold,
            thread_count=cfg["num_threads"],
            verbose=False,
            allow_writing_files=False,
        ),
        "member_extratrees": ExtraTreesRegressor(
            n_estimators=500,
            min_samples_leaf=12,
            max_features=0.8,
            n_jobs=cfg["num_threads"],
            random_state=42 + fold,
        ),
    }
    for method, model in models.items():
        scalar = method.startswith("scalar_")
        x = x_scalar if scalar else x_members
        if method == "scalar_isotonic":
            model.fit(x[train, 0], y[train])
            prediction = model.predict(x[valid, 0])
        else:
            model.fit(x[train], y[train])
            prediction = model.predict(x[valid])
        assert np.isfinite(prediction).all()
        rows[method].append(
            data.filter(pl.Series(valid))
            .select(
                "compound_id",
                "endpoint",
                "fold",
                "is_neighbor",
                "conf_low",
                "conf_high",
                "y_true",
            )
            .with_columns(
                pl.lit(method).alias("method"),
                pl.lit("tdi_probability_crossfit").alias("split_id"),
                pl.Series("prediction", prediction),
            )
        )
        fits.append(
            {"method": method, "fold": fold, "train": int(train.sum()), "valid": int(valid.sum())}
        )

metrics = []
for method, frames in rows.items():
    frame = pl.concat(frames).sort("fold", "compound_id")
    path = out / f"{method}-oof.parquet"
    frame.write_parquet(path)
    score = float(
        scorer.rae_soft_threshold_absolute_error(
            frame["y_true"].to_numpy(),
            frame["prediction"].to_numpy(),
            y_true_lower=frame["conf_low"].to_numpy(),
            y_true_upper=frame["conf_high"].to_numpy(),
        )
    )
    metrics.append(
        {
            "method": method,
            "ST_RAE": score,
            "MAE": float(np.abs(frame["y_true"] - frame["prediction"]).mean()),
            "Spearman": float(
                frame.select(pl.corr("y_true", "prediction", method="spearman")).item()
            ),
        }
    )

metrics = sorted(metrics, key=lambda item: item["ST_RAE"])
pl.DataFrame(metrics).write_csv(out / "metrics.csv")
pl.DataFrame(fits).write_csv(out / "fits.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "config": cfg,
            "rows": data.height,
            "member_columns": member_columns,
            "tdi_pIC50_pearson": float(data.select(pl.corr("y_true", "tdi_equal7")).item()),
            "tdi_pIC50_spearman": float(
                data.select(pl.corr("y_true", "tdi_equal7", method="spearman")).item()
            ),
            "metrics": metrics,
        },
        indent=2,
    )
    + "\n"
)
print(json.dumps(metrics), flush=True)
