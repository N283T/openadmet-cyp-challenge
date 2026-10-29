"""CYP3A4 direct-label augmentation with TDI-condition pIC50 of TDI-negative compounds.

The challenge TDI file carries a TDI-condition 3A4 pIC50 for 1,249 compounds that
have no direct-inhibition label. For TDI-negative compounds the two conditions
differ by less than two-fold by definition, so a shifted TDI-condition value is a
proxy direct label. Each recipe is fitted twice on identical features and folds,
without and with proxy rows, so the paired difference isolates the augmentation.

Proxy rows enter a fold's training set only when their compound's parent fold is a
training fold; the proxy shift is re-estimated from that fold's training compounds
that carry both labels. Validation rows and folds are the incumbent's OOF rows.
"""

import importlib.util
import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
import torch
from tabpfn import TabPFNRegressor
from tabpfn.constants import ModelVersion
from threadpoolctl import threadpool_limits

DIRECT = "CYP3A4_pIC50_direct_inhibition"
TDI = "CYP3A4_pIC50_TDI_condition"
LABEL = "CYP3A4_is_TDI"
AUX = [f"aux_{i}" for i in range(4)]


def matrix(ids, features, columns, nested_fold):
    """Feature matrix keyed by compound ID, with optional fold-specific nested log2fc."""
    table = ids.join(features, on="compound_id", how="left", maintain_order="left")
    x = table.select(columns).to_numpy().astype(np.float32)
    if nested_fold is not None:
        extra = ids.join(nested_fold, on="compound_id", how="left", maintain_order="left")
        x = np.column_stack([x, extra.select(AUX).to_numpy()])
    assert np.isfinite(x).all()
    return x


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    threadpool_limits(limits=cfg["threads"])
    torch.set_num_threads(cfg["threads"])
    spec = importlib.util.spec_from_file_location("official", inp / cfg["scorer"])
    scorer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scorer)

    def st_rae(frame, prediction):
        return float(
            scorer.rae_soft_threshold_absolute_error(
                frame["y_true"].to_numpy(),
                prediction,
                y_true_upper=frame["conf_high"].to_numpy(),
                y_true_lower=frame["conf_low"].to_numpy(),
            )
        )

    measurements = pl.read_parquet(inp / "measurements.parquet")
    wide = measurements.filter(pl.col("endpoint").is_in([DIRECT, TDI, LABEL])).pivot(
        on="endpoint", index="compound_id", values="value", aggregate_function="first"
    )
    direct = measurements.filter(pl.col("endpoint") == DIRECT).select(
        "compound_id", "value", "conf_low", "conf_high"
    )
    assignments = pl.read_parquet(inp / "parent-assignments.parquet").filter(
        pl.col("strategy") == "series"
    )
    control = (
        pl.read_parquet(inp / "control-3a4.parquet")
        .filter((pl.col("method") == "caruana") & (pl.col("endpoint") == DIRECT))
        .sort("compound_id")
    )
    assert control.height == cfg["expected_validation_rows"]

    proxy = (
        wide.filter(pl.col(DIRECT).is_null() & pl.col(TDI).is_not_null() & (pl.col(LABEL) == 0))
        .join(assignments.select("compound_id", "fold"), on="compound_id", validate="1:1")
        .select("compound_id", "fold", pl.col(TDI).alias("tdi"))
    )
    paired = (
        wide.filter(pl.col(DIRECT).is_not_null() & pl.col(TDI).is_not_null() & (pl.col(LABEL) == 0))
        .join(assignments.select("compound_id", "fold"), on="compound_id", validate="1:1")
        .select("compound_id", "fold", (pl.col(TDI) - pl.col(DIRECT)).alias("shift"))
    )
    assert proxy.height == cfg["expected_proxy_rows"]
    train_pool = direct.join(
        assignments.select("compound_id", "fold", "is_anchor"), on="compound_id", validate="1:1"
    )

    def model(readout, fold):
        if readout == "tabpfn35":
            return TabPFNRegressor.create_default_for_version(
                ModelVersion(cfg["tabpfn_version"]),
                n_estimators=cfg["tabpfn_estimators"],
                model_path=str(inp / "tabpfn.safetensors"),
                device="cuda",
                random_state=42 + fold,
                auto_scale_n_estimators=False,
                n_preprocessing_jobs=1,
                memory_saving_mode="auto",
                show_progress_bar=False,
                ignore_pretraining_limits=True,
                fit_mode="fit_preprocessors",
            )
        return lgb.LGBMRegressor(
            objective="regression_l1",
            n_estimators=300,
            num_leaves=15,
            min_child_samples=20,
            learning_rate=0.05,
            reg_lambda=1.0,
            n_jobs=cfg["threads"],
            random_state=42 + fold,
            verbosity=-1,
            deterministic=True,
            force_col_wise=True,
        )

    nested = pl.read_parquet(inp / "nested-log2fc.parquet")
    incumbent = control["prediction"].to_numpy()
    records, oof_frames, shifts = [], [], []
    for recipe in cfg["recipes"]:
        features = pl.read_parquet(inp / recipe["file"])
        columns = sorted(c for c in features.columns if c.startswith(tuple(recipe["prefixes"])))
        assert columns and features["compound_id"].n_unique() == features.height
        predictions = {}
        for arm in ("baseline", "proxy"):
            member = np.full(control.height, np.nan)
            for fold in range(5):
                train = train_pool.filter(pl.col("fold") != fold).select(
                    "compound_id", pl.col("value").alias("y"), pl.lit(1.0).alias("w")
                )
                if arm == "proxy":
                    shift = float(paired.filter(pl.col("fold") != fold)["shift"].median())
                    shifts.append({"recipe": recipe["name"], "fold": fold, "shift": shift})
                    train = pl.concat(
                        [
                            train,
                            proxy.filter(pl.col("fold") != fold).select(
                                "compound_id",
                                (pl.col("tdi") - shift).alias("y"),
                                pl.lit(cfg["proxy_weight"]).alias("w"),
                            ),
                        ]
                    )
                valid = control.filter(pl.col("fold") == fold)

                nested_fold = (
                    nested.filter(pl.col("outer_fold") == fold).select("compound_id", *AUX)
                    if recipe.get("nested")
                    else None
                )
                fitted = model(recipe["readout"], fold)
                xt = matrix(train.select("compound_id"), features, columns, nested_fold)
                if recipe["readout"] == "lgbm":
                    fitted.fit(xt, train["y"].to_numpy(), sample_weight=train["w"].to_numpy())
                else:
                    fitted.fit(xt, train["y"].to_numpy())
                mask = control["fold"].to_numpy() == fold
                member[mask] = fitted.predict(
                    matrix(valid.select("compound_id"), features, columns, nested_fold)
                )
                del fitted
                torch.cuda.empty_cache()
            assert np.isfinite(member).all()
            predictions[arm] = member
            oof_frames.append(
                control.select(
                    "compound_id", "endpoint", "fold", "y_true", "conf_low", "conf_high"
                ).with_columns(
                    pl.lit(f"{recipe['name']}--{arm}").alias("method"),
                    pl.lit(cfg["split_id"]).alias("split_id"),
                    pl.Series("prediction", member),
                )
            )
        folds = control["fold"].to_numpy()
        record = {"recipe": recipe["name"], "incumbent": st_rae(control, incumbent)}
        for arm, member in predictions.items():
            record[f"{arm}_st_rae"] = st_rae(control, member)
            for weight in cfg["blend_weights"]:
                record[f"{arm}_blend_{weight}"] = st_rae(
                    control, (1 - weight) * incumbent + weight * member
                )
        record["proxy_fold_wins"] = int(
            sum(
                st_rae(control.filter(pl.Series(folds == f)), predictions["proxy"][folds == f])
                < st_rae(control.filter(pl.Series(folds == f)), predictions["baseline"][folds == f])
                for f in range(5)
            )
        )
        records.append(record)
        print(json.dumps(record), flush=True)
        pl.DataFrame(records).write_parquet(out / "paired-metrics.parquet")
        pl.concat(oof_frames).write_parquet(out / "oof.parquet")
    pl.DataFrame(shifts).write_parquet(out / "proxy-shifts.parquet")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "proxy_rows": proxy.height,
                "paired_rows": paired.height,
                "proxy_weight": cfg["proxy_weight"],
                "records": records,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
