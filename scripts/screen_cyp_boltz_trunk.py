"""Fixed five-fold PXR-style trunk ablations; exploratory, no candidate selection."""

from __future__ import annotations

import gc
import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
import torch
from lightgbm import LGBMClassifier, LGBMRegressor
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score
from tabpfn import TabPFNClassifier, TabPFNRegressor
from tabpfn.constants import ModelVersion
from threadpoolctl import threadpool_limits


def feature_views(columns):
    def choose(block, stats):
        return [c for c in columns if any(c.startswith(f"{block}_{stat}_") for stat in stats)]

    base = choose("s_protein", ["mean"]) + choose("s_ligand", ["mean"])
    pair = lambda region, direction="pl", stats=("mean", "max"): choose(  # noqa: E731
        f"z_{region}_{direction}", stats
    )
    views = {"ligand_only": choose("s_ligand", ["mean"]) + choose("z_ll", ["mean", "max"])}
    for region in ["protein", "core", "extended"]:
        views[f"pxr_{region}"] = base + pair(region)
    for region in ["core", "extended"]:
        views[f"{region}_mixed_hem"] = base + pair(f"{region}_hem")
        views[f"{region}_separate_hem"] = base + pair(region) + pair("hem")
    views["core_separate_fe"] = base + pair("core") + pair("hem_other") + pair("fe")
    views["core_bidirectional_hem"] = (
        base + pair("core") + pair("core", "lp") + pair("hem") + pair("hem", "lp")
    )
    views["core_separate_hem_s"] = views["core_separate_hem"] + choose("s_hem", ["mean"])
    views["core_distribution_hem"] = (
        base
        + choose("s_ligand", ["std"])
        + pair("core", stats=("mean", "std", "q10", "q90"))
        + pair("hem", stats=("mean", "std", "q10", "q90"))
    )
    return views


def metrics(frame, p, classification):
    y = frame["value"].to_numpy()
    if classification:
        return {
            "AP": float(average_precision_score(y, p)),
            "AUROC": float(roc_auc_score(y, p)),
            "MCC_at_0.5": float(matthews_corrcoef(y, p >= 0.5)),
        }
    low, high = frame["conf_low"].to_numpy(), frame["conf_high"].to_numpy()
    error = lambda prediction: np.maximum(low - prediction, 0) + np.maximum(prediction - high, 0)  # noqa: E731
    return {
        "ST_RAE": float(error(p).sum() / error(y.mean()).sum()),
        "MAE": float(np.abs(y - p).mean()),
        "Spearman": float(spearmanr(y, p).statistic),
    }


def tree(classification, method, seed, threads):
    cls = LGBMClassifier if classification else LGBMRegressor
    return cls(
        objective="binary"
        if classification
        else ("regression_l1" if method == "lgbm_l1" else "regression"),
        n_estimators=400,
        num_leaves=15,
        min_child_samples=20,
        learning_rate=0.04,
        reg_lambda=5,
        colsample_bytree=0.8,
        n_jobs=threads,
        verbosity=-1,
        deterministic=True,
        force_col_wise=True,
        random_state=seed,
    )


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    threadpool_limits(limits=cfg["threads"])
    torch.set_num_threads(cfg["threads"])
    classification = cfg["track"] == "classification"
    labels = (
        pl.read_parquet(inp / "measurements.parquet")
        .filter(pl.col("endpoint") == cfg["endpoint"])
        .select("compound_id", "value", "conf_low", "conf_high")
    )
    if cfg.get("cohort_from_control"):
        labels = pl.read_parquet(inp / "control-oof.parquet")
        for column, value in cfg["control_filter"].items():
            labels = labels.filter(pl.col(column) == value)
        labels = labels.filter(pl.col("endpoint") == cfg["endpoint"])
        if labels.height != cfg["expected_labelled_rows"]:
            raise ValueError("Incumbent cohort count differs from the declared protocol")
        labels = labels.select(
            "compound_id",
            pl.col("y_true").alias("value"),
            "conf_low",
            "conf_high",
            pl.col("fold").alias("control_fold"),
        )
    folds = pl.read_parquet(inp / "assignments.parquet").filter(
        pl.col("strategy") == cfg.get("assignment_strategy", "umap")
    )
    features = pl.read_parquet(inp / "trunk.parquet")
    cheme = pl.concat([pl.read_parquet(p) for p in sorted(inp.glob("cheme-*.parquet"))])
    cheme_columns = [c for c in cheme.columns if c.startswith("chemeleon_")]
    mordred = pl.concat([pl.read_parquet(p) for p in sorted(inp.glob("mordred-*.parquet"))])
    mordred_columns = [c for c in mordred.columns if c != "compound_id"]
    cheme = cheme.join(mordred, on="compound_id", validate="1:1")
    if len(cheme_columns) != 2048:
        raise ValueError("Incomplete CheMeleon columns")
    cheme_columns += mordred_columns
    views = feature_views(features.columns)
    if set(views) != set(cfg["views"]):
        raise ValueError("Configured views differ from implementation")
    (out / "feature_views.json").write_text(json.dumps(views, indent=2) + "\n")
    frame = (
        features.join(cheme, on="compound_id", validate="1:1")
        .join(
            labels,
            on="compound_id",
            how="left",
            validate="1:1",
        )
        .join(
            folds.select("compound_id", "fold", "group"),
            on="compound_id",
            how="left",
            validate="1:1",
        )
        .sort("compound_id")
    )
    labelled = frame["value"].is_not_null().to_numpy() & (frame["role"].to_numpy() == "train")
    if cfg.get("cohort_from_control"):
        matched = frame.filter(pl.Series(labelled))
        if not matched["fold"].equals(matched["control_fold"]):
            raise ValueError("Trunk folds differ from incumbent OOF folds")
        matched.select(
            "compound_id", "fold", "group", "value", "conf_low", "conf_high"
        ).write_parquet(out / "evaluation-cohort.parquet")
    y = frame["value"].to_numpy()
    f = frame["fold"].fill_null(-1).to_numpy()
    test = frame["role"].to_numpy() == "test"
    if test.sum() != 750 or (labelled & (f < 0)).any():
        raise ValueError("Missing test rows or labelled fold assignments")
    expected = labels["compound_id"].to_list()
    missing = sorted(set(expected) - set(frame.filter(pl.Series(labelled))["compound_id"]))
    if missing:
        raise ValueError(f"Missing labelled molecules: {len(missing)}")
    recipes = [("molecular_only", "molecular", cheme_columns)]
    for view, cols in views.items():
        recipes += [(view, "trunk", cols), (view, "trunk_molecular", cols + cheme_columns)]
    methods = (
        ["lgbm_binary", "tabpfn_top500"]
        if classification
        else ["lgbm_l1", "lgbm_l2", "tabpfn_top500"]
    )
    (out / "models").mkdir()
    (out / "predictions").mkdir()
    rows, aggregate, predictions, failures = [], [], [], []
    baseline_predictions = {}
    for view, mode, cols in recipes:
        x = frame.select(cols).to_numpy().astype(np.float32)
        x[~np.isfinite(x)] = np.nan
        for method in methods:
            name = f"{view}--{mode}--{method}"
            oof = np.full(len(frame), np.nan)
            test_predictions = []
            for fold in range(5):
                train, valid = labelled & (f != fold), labelled & (f == fold)
                held_pool = (f == fold) & ~test
                predict_mask = held_pool | test
                seed = cfg["seed"] + fold
                train_groups = set(frame.filter(pl.Series(train))["group"])
                valid_groups = set(frame.filter(pl.Series(held_pool))["group"])
                if not train_groups.isdisjoint(valid_groups):
                    raise ValueError("Group crosses the training/held-out boundary")
                model_path = out / "models" / f"{name}--fold{fold}"
                try:
                    medians = np.nanmedian(x[train], axis=0)
                    medians[~np.isfinite(medians)] = 0.0
                    xx = np.where(np.isfinite(x), x, medians).astype(np.float32)
                    selected = np.arange(len(cols))
                    selector = None
                    if method == "tabpfn_top500" and len(cols) > 500:
                        selector = tree(classification, "lgbm_l1", seed, cfg["threads"])
                        selector.fit(xx[train], y[train])
                        gain = selector.booster_.feature_importance(importance_type="gain")
                        selected = np.lexsort((np.array(cols), -gain))[:500]
                    if method == "tabpfn_top500":
                        cls = TabPFNClassifier if classification else TabPFNRegressor
                        model = cls.create_default_for_version(
                            ModelVersion.V2_6,
                            n_estimators=cfg["tabpfn_estimators"],
                            auto_scale_n_estimators=False,
                            model_path=str(inp / "tabpfn.ckpt"),
                            device="cuda",
                            random_state=seed,
                            n_preprocessing_jobs=1,
                            memory_saving_mode="auto",
                            show_progress_bar=False,
                            ignore_pretraining_limits=True,
                            fit_mode="fit_preprocessors",
                        )
                    else:
                        model = tree(classification, method, seed, cfg["threads"])
                    model.fit(xx[train][:, selected], y[train])
                    Path(str(model_path) + ".parameters.json").write_text(
                        json.dumps(model.get_params(), indent=2, default=str) + "\n"
                    )
                    query = xx[predict_mask][:, selected]
                    p = model.predict_proba(query)[:, 1] if classification else model.predict(query)
                    p = np.asarray(p, dtype=np.float64)
                    if not np.isfinite(p).all():
                        raise ValueError("Nonfinite predictions")
                    indices = np.flatnonzero(predict_mask)
                    full = np.full(len(frame), np.nan)
                    full[indices] = p
                    oof[held_pool] = full[held_pool]
                    test_predictions.append(full[test])
                    pred = (
                        frame.filter(pl.Series(predict_mask))
                        .select("compound_id", "group", "role")
                        .with_columns(
                            pl.lit(cfg["endpoint"]).alias("endpoint"),
                            pl.lit(cfg["split_id"]).alias("split_id"),
                            pl.lit(fold).alias("fold"),
                            pl.lit(name).alias("model"),
                            pl.Series("prediction", p),
                        )
                    )
                    pred.write_parquet(out / "predictions" / f"{name}--fold{fold}.parquet")
                    predictions.append(pred)
                    row = {
                        "view": view,
                        "mode": mode,
                        "method": method,
                        "fold": fold,
                        "n_train": int(train.sum()),
                        "n_valid": int(valid.sum()),
                        "n_features": len(selected),
                        **metrics(frame.filter(pl.Series(valid)), full[valid], classification),
                    }
                    rows.append(row)
                    pl.DataFrame(rows).write_parquet(out / "fold-metrics.parquet")
                    joblib.dump(
                        {
                            "columns": [cols[i] for i in selected],
                            "selector": selector,
                            "imputation_medians": medians[selected],
                            "train_ids": frame.filter(pl.Series(train))["compound_id"].to_list(),
                            "fold": fold,
                            "seed": seed,
                        },
                        str(model_path) + ".preprocessing.joblib",
                    )
                    if method == "tabpfn_top500":
                        model.save_fit_state(str(model_path) + ".tabpfn_fit")
                    else:
                        joblib.dump(model, str(model_path) + ".joblib")
                    del model, selector
                    gc.collect()
                    torch.cuda.empty_cache()
                    print(json.dumps({"model": name, **row}), flush=True)
                except Exception as exc:
                    failure = {"model": name, "fold": fold, "error": repr(exc)}
                    failures.append(failure)
                    (out / "failures.json").write_text(json.dumps(failures, indent=2) + "\n")
                    print(json.dumps(failure), flush=True)
                    gc.collect()
                    torch.cuda.empty_cache()
                    break
            if not np.isfinite(oof[labelled]).all() or len(test_predictions) != 5:
                continue
            record = {
                "view": view,
                "mode": mode,
                "method": method,
                "model": name,
                **metrics(frame.filter(pl.Series(labelled)), oof[labelled], classification),
            }
            if mode == "molecular":
                baseline_predictions[method] = oof
            control = baseline_predictions[method]
            if not classification:
                record["residual_correlation_vs_molecular"] = float(
                    np.corrcoef(y[labelled] - oof[labelled], y[labelled] - control[labelled])[0, 1]
                )
                for weight in [0.1, 0.25, 0.5]:
                    record[f"blend_{weight}_ST_RAE"] = metrics(
                        frame.filter(pl.Series(labelled)),
                        (1 - weight) * control[labelled] + weight * oof[labelled],
                        False,
                    )["ST_RAE"]
            aggregate.append(record)
            pl.DataFrame(aggregate).write_parquet(out / "aggregate-metrics.parquet")
            test_frame = (
                frame.filter(pl.Series(test))
                .select("compound_id")
                .with_columns(
                    pl.lit(cfg["endpoint"]).alias("endpoint"),
                    pl.lit(cfg["split_id"]).alias("split_id"),
                    pl.lit(-1).alias("fold"),
                    pl.lit(name).alias("model"),
                    pl.Series("prediction", np.mean(test_predictions, axis=0)),
                )
            )
            test_frame.write_parquet(out / "predictions" / f"{name}--test-mean.parquet")
    if predictions:
        pl.concat(predictions).write_parquet(out / "predictions.parquet")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "endpoint": cfg["endpoint"],
                "successful_recipes": len(aggregate),
                "planned_recipes": len(recipes) * len(methods),
                "failures": failures,
                "exploratory_only": True,
                "classification_deployment_supported": False,
                "classification_limit": (
                    "Target-extracted pool does not cover all 6146 molecules; "
                    "no call-rate selection"
                ),
                "blends": (
                    "Prespecified weights vs matched molecular-only control; "
                    "exploratory, no selection"
                ),
                "pool_train_rows": int((~test).sum()),
                "labelled_rows": int(labelled.sum()),
            },
            indent=2,
        )
        + "\n"
    )
    if failures:
        raise RuntimeError(f"{len(failures)} recipes failed; successful outputs preserved")


if __name__ == "__main__":
    main()
