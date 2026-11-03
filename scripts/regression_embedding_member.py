"""Frozen-embedding regression members on the parent cohort, scored against the incumbent.

Each endpoint uses the incumbent ensemble's own OOF rows and folds, so member OOF,
incumbent OOF and fixed-weight blends are directly comparable. Every fold model also
predicts the 750 blind compounds; the blind member is the five-fold mean.
"""

import gc
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
from lightgbm import LGBMRegressor
from scipy.stats import spearmanr
from tabicl import TabICLRegressor
from tabpfn import TabPFNRegressor
from tabpfn.constants import ModelVersion
from threadpoolctl import threadpool_limits


def metrics(y, low, high, p):
    def error(q):
        return np.maximum(low - q, 0) + np.maximum(q - high, 0)

    return {
        "ST_RAE": float(error(p).sum() / error(y.mean()).sum()),
        "MAE": float(np.abs(y - p).mean()),
        "Spearman": float(spearmanr(y, p).statistic),
    }


def lgbm(seed, threads):
    return LGBMRegressor(
        objective="regression_l1",
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


def build_model(readout, seed, cfg, inp):
    if readout == "tabpfn":
        return TabPFNRegressor.create_default_for_version(
            ModelVersion(cfg.get("tabpfn_version", "v2.6")),
            n_estimators=cfg["tabpfn_estimators"],
            auto_scale_n_estimators=False,
            model_path=str(inp / cfg.get("tabpfn_weights", "tabpfn.ckpt")),
            device="cuda",
            random_state=seed,
            n_preprocessing_jobs=1,
            memory_saving_mode="auto",
            show_progress_bar=False,
            ignore_pretraining_limits=True,
            fit_mode="fit_preprocessors",
        )
    if readout == "tabicl":
        return TabICLRegressor(
            n_estimators=cfg["tabicl_estimators"],
            batch_size=1,
            device="cuda",
            n_jobs=cfg["threads"],
            model_path=str(inp / "tabicl.ckpt"),
            allow_auto_download=False,
            use_amp="auto",
            offload_mode="auto",
            disk_offload_dir=os.environ["TMPDIR"],
            random_state=seed,
        )
    if readout == "lgbm_l1":
        return lgbm(seed, cfg["threads"])
    raise ValueError(readout)


def select_columns(x_train, y_train, selection, seed, threads):
    if selection == "all":
        return np.arange(x_train.shape[1])
    kind, k = selection.split(":")
    k = min(int(k), x_train.shape[1])
    if kind == "variance":
        return np.sort(np.argsort(-x_train.var(axis=0), kind="stable")[:k])
    if kind == "gain":
        selector = lgbm(seed, threads).fit(x_train, y_train)
        gain = selector.booster_.feature_importance(importance_type="gain")
        return np.sort(np.argsort(-gain, kind="stable")[:k])
    raise ValueError(selection)


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    threadpool_limits(limits=cfg["threads"])
    torch.set_num_threads(cfg["threads"])
    (out / "predictions").mkdir()

    blocks = {}
    for name, block in cfg["feature_blocks"].items():
        frame = pl.concat([pl.read_parquet(inp / f) for f in block["files"]], how="vertical")
        cols = [c for c in frame.columns if c.startswith(tuple(block["prefixes"]))]
        assert cols and frame["compound_id"].n_unique() == frame.height
        blocks[name] = frame.select("compound_id", *cols).rename({c: f"{name}:{c}" for c in cols})
    blind_ids = pl.read_parquet(inp / cfg["blind_compounds"])["compound_id"].sort()
    assert blind_ids.n_unique() == blind_ids.len() == 750

    aggregate, fold_rows, oof_frames, blind_frames = [], [], [], []
    for endpoint, control in cfg["controls"].items():
        oof = pl.read_parquet(inp / control["file"])
        for column, value in control["filter"].items():
            oof = oof.filter(pl.col(column) == value)
        oof = oof.filter(pl.col("endpoint") == endpoint).sort("compound_id")
        assert oof.height == control["rows"] and oof["compound_id"].n_unique() == oof.height
        y, low, high = (oof[c].to_numpy() for c in ("y_true", "conf_low", "conf_high"))
        folds = oof["fold"].to_numpy()
        incumbent = oof["prediction"].to_numpy()
        reference = metrics(y, low, high, incumbent)
        for recipe in cfg["recipes"]:
            name = f"{recipe['blocks']}--{recipe['selection']}--{recipe['readout']}"
            ids = pl.concat([oof.select("compound_id"), pl.DataFrame({"compound_id": blind_ids})])
            table = ids
            for block in recipe["blocks"].split("+"):
                table = table.join(blocks[block], on="compound_id", how="left", validate="1:1")
            assert table["compound_id"].equals(ids["compound_id"])
            x = table.drop("compound_id").to_numpy().astype(np.float32)
            x[~np.isfinite(x)] = np.nan
            n = oof.height
            member = np.full(n, np.nan)
            blind = []
            for fold in range(5):
                train, valid = folds != fold, folds == fold
                seed = cfg["seed"] + fold
                medians = np.nanmedian(x[:n][train], axis=0)
                medians[~np.isfinite(medians)] = 0.0
                xx = np.where(np.isfinite(x), x, medians).astype(np.float32)
                columns = select_columns(
                    xx[:n][train], y[train], recipe["selection"], seed, cfg["threads"]
                )
                model = build_model(recipe["readout"], seed, cfg, inp)
                model.fit(xx[:n][train][:, columns], y[train])
                member[valid] = np.asarray(model.predict(xx[:n][valid][:, columns]), dtype=float)
                blind.append(np.asarray(model.predict(xx[n:][:, columns]), dtype=float))
                fold_rows.append(
                    {
                        "endpoint": endpoint,
                        "model": name,
                        "fold": fold,
                        "n_train": int(train.sum()),
                        "n_features": len(columns),
                        **metrics(y[valid], low[valid], high[valid], member[valid]),
                    }
                )
                print(json.dumps(fold_rows[-1]), flush=True)
                del model
                gc.collect()
                torch.cuda.empty_cache()
            assert np.isfinite(member).all() and np.isfinite(blind).all()
            record = {
                "endpoint": endpoint,
                "model": name,
                **recipe,
                **metrics(y, low, high, member),
                "incumbent_ST_RAE": reference["ST_RAE"],
                "member_incumbent_corr": float(np.corrcoef(member, incumbent)[0, 1]),
                "residual_corr": float(np.corrcoef(y - member, y - incumbent)[0, 1]),
            }
            for weight in cfg["blend_weights"]:
                blended = (1 - weight) * incumbent + weight * member
                record[f"blend_{weight}_ST_RAE"] = metrics(y, low, high, blended)["ST_RAE"]
                record[f"blend_{weight}_fold_wins"] = int(
                    sum(
                        metrics(y[f], low[f], high[f], blended[f])["ST_RAE"]
                        < metrics(y[f], low[f], high[f], incumbent[f])["ST_RAE"]
                        for f in (folds == k for k in range(5))
                    )
                )
            aggregate.append(record)
            print(json.dumps(record), flush=True)
            oof_frames.append(
                oof.select(
                    "compound_id", "endpoint", "fold", "y_true", "conf_low", "conf_high"
                ).with_columns(
                    pl.lit(cfg["split_id"]).alias("split_id"),
                    pl.lit(name).alias("model"),
                    pl.Series("prediction", member),
                )
            )
            blind_frames.append(
                pl.DataFrame({"compound_id": blind_ids}).with_columns(
                    pl.lit(endpoint).alias("endpoint"),
                    pl.lit(name).alias("model"),
                    pl.Series("prediction", np.mean(blind, axis=0)),
                )
            )
            pl.DataFrame(aggregate).write_parquet(out / "aggregate-metrics.parquet")
            pl.DataFrame(fold_rows).write_parquet(out / "fold-metrics.parquet")
            pl.concat(oof_frames).write_parquet(out / "predictions/oof.parquet")
            pl.concat(blind_frames).write_parquet(out / "predictions/blind-mean.parquet")
    (out / "summary.json").write_text(
        json.dumps(
            {
                "recipes": len(cfg["recipes"]),
                "endpoints": list(cfg["controls"]),
                "completed": len(aggregate),
                "selection": "Exploratory grid; blends use prespecified fixed weights",
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
