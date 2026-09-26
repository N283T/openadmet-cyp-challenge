"""Compare D6 potency mapping with measured and fold-excluded primary response."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import spearmanr
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LinearRegression


def st_rae(y: np.ndarray, lo: np.ndarray, hi: np.ndarray, pred: np.ndarray) -> float:
    center = y.mean()
    numerator = (np.maximum(lo - pred, 0) + np.maximum(pred - hi, 0)).sum()
    denominator = (np.maximum(lo - center, 0) + np.maximum(center - hi, 0)).sum()
    return float(numerator / denominator)


def mapped(name: str, train_x: np.ndarray, train_y: np.ndarray, test_x: np.ndarray) -> np.ndarray:
    if name == "linear":
        model = LinearRegression().fit(train_x.reshape(-1, 1), train_y)
        return model.predict(test_x.reshape(-1, 1))
    if name == "isotonic":
        model = IsotonicRegression(increasing=False, out_of_bounds="clip")
        return model.fit(train_x, train_y).predict(test_x)
    raise ValueError(name)


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    direct = (
        pl.read_parquet(inp / "measurements.parquet")
        .filter(pl.col("endpoint") == cfg["direct_endpoint"])
        .select("compound_id", pl.col("value").alias("y_true"), "conf_low", "conf_high")
    )
    primary = (
        pl.read_parquet(inp / "regression-labels.parquet")
        .filter(pl.col("head") == cfg["primary_head"])
        .select("compound_id", pl.col("value").alias("measured_primary"))
    )
    assignment = pl.read_parquet(inp / "parent-assignments.parquet").select(
        "compound_id", "fold", "is_anchor"
    )
    nested = pl.read_parquet(inp / "nested-primary.parquet").select(
        "compound_id", "outer_fold", pl.col(cfg["primary_feature"]).alias("predicted_primary")
    )
    control = (
        pl.read_parquet(inp / "control-oof.parquet")
        .filter(pl.col("endpoint") == cfg["direct_endpoint"])
        .select("compound_id", "fold", "prediction")
    )
    labels = assignment.join(direct, on="compound_id", how="inner", validate="1:1").join(
        primary, on="compound_id", validate="1:1"
    )
    if labels.height != 1493 or control.height != 1457:
        raise ValueError("Unexpected D6 label or control coverage")
    results = []
    associations = []
    for fold in range(5):
        frame = labels.join(
            nested.filter(pl.col("outer_fold") == fold).drop("outer_fold"),
            on="compound_id",
            validate="1:1",
        ).sort("compound_id")
        train = frame.filter((pl.col("fold") != fold) | pl.col("is_anchor"))
        valid = frame.filter((pl.col("fold") == fold) & ~pl.col("is_anchor"))
        if valid.height < 100 or train.height + valid.height != 1493:
            raise ValueError(f"Incomplete fold {fold}")
        fit_y = train["y_true"].to_numpy()
        truth = valid["y_true"].to_numpy()
        measured = valid["measured_primary"].to_numpy()
        predicted = valid["predicted_primary"].to_numpy()
        associations.append(
            {
                "fold": fold,
                "n": valid.height,
                "measured_primary_spearman": float(spearmanr(truth, measured).statistic),
                "predicted_primary_spearman": float(spearmanr(truth, predicted).statistic),
                "primary_prediction_rmse": float(np.mean((predicted - measured) ** 2) ** 0.5),
            }
        )
        for feature, train_col, valid_col in [
            ("measured_oracle", "measured_primary", "measured_primary"),
            ("predicted", "predicted_primary", "predicted_primary"),
        ]:
            for mapping in cfg["mappings"]:
                prediction = mapped(
                    mapping, train[train_col].to_numpy(), fit_y, valid[valid_col].to_numpy()
                )
                results.append(
                    valid.select(
                        "compound_id", "fold", "y_true", "conf_low", "conf_high"
                    ).with_columns(
                        pl.lit(cfg["direct_endpoint"]).alias("endpoint"),
                        pl.lit(cfg["split"]).alias("split_id"),
                        pl.lit(f"{feature}_{mapping}").alias("method"),
                        pl.Series("prediction", prediction),
                    )
                )
        control_part = valid.select("compound_id", "fold", "y_true", "conf_low", "conf_high").join(
            control.filter(pl.col("fold") == fold), on=["compound_id", "fold"], validate="1:1"
        )
        if control_part.height != valid.height:
            raise ValueError(f"Incomplete control fold {fold}")
        results.append(
            control_part.with_columns(
                pl.lit(cfg["direct_endpoint"]).alias("endpoint"),
                pl.lit(cfg["split"]).alias("split_id"),
                pl.lit("submitted_control").alias("method"),
            )
        )
    columns = [
        "compound_id",
        "fold",
        "y_true",
        "conf_low",
        "conf_high",
        "endpoint",
        "split_id",
        "method",
        "prediction",
    ]
    oof = pl.concat([frame.select(columns) for frame in results]).sort("method", "compound_id")
    if oof.height != 5 * 1457 or oof.select("method", "compound_id").is_duplicated().any():
        raise ValueError("Incomplete diagnostic OOF")
    metrics = []
    for method, frame in oof.group_by("method"):
        for fold in [-1, 0, 1, 2, 3, 4]:
            part = frame if fold == -1 else frame.filter(pl.col("fold") == fold)
            y, lo, hi, pred = (
                part[c].to_numpy() for c in ["y_true", "conf_low", "conf_high", "prediction"]
            )
            metrics.append(
                {
                    "method": method[0],
                    "fold": fold,
                    "n": part.height,
                    "ST_RAE": st_rae(y, lo, hi, pred),
                    "MAE": float(np.abs(y - pred).mean()),
                    "Spearman": float(spearmanr(y, pred).statistic),
                }
            )
    oof.write_parquet(out / "oof.parquet")
    pl.DataFrame(metrics).sort("method", "fold").write_csv(out / "metrics.csv")
    pl.DataFrame(associations).write_csv(out / "associations.csv")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "config": cfg,
                "oof_rows": oof.height,
                "qualification": "Measured-primary arm uses unavailable test-time measurements",
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
