"""Compare Mordred alone and appended to the saved fold-safe 260-d features."""

import gc
import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
import torch
from sklearn.impute import SimpleImputer
from split_baselines import metrics
from tabpfn import TabPFNClassifier, TabPFNRegressor
from tabpfn.constants import ModelVersion

from cyp_challenge.storage import write_json


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    torch.set_num_threads(4)
    torch.manual_seed(cfg["seed"])
    fold = cfg["fold"]
    labels = pl.read_parquet(next(inp.glob("*-measurements.parquet")))
    folds = pl.read_parquet(next(inp.glob("*-assignments.parquet"))).filter(
        pl.col("strategy") == "umap"
    )
    desc = pl.concat([pl.read_parquet(p) for p in sorted(inp.glob("*-mordred-*.parquet"))])
    columns = [c for c in desc.columns if c != "compound_id"]
    all_rows, scores = [], []
    for ep in cfg["endpoints"]:
        frame = (
            labels.filter(pl.col("endpoint") == ep)
            .join(folds.select("compound_id", "fold"), on="compound_id", validate="1:1")
            .join(desc, on="compound_id", how="left", validate="1:1")
            .sort("compound_id")
        )
        train = frame["fold"].to_numpy() != fold
        y = frame["value"].to_numpy()
        x = frame.select(columns).to_numpy().astype(np.float64)
        x[~np.isfinite(x)] = np.nan
        keep = np.isfinite(x[train]).any(axis=0)
        # Drop only train-empty/constant columns; no label-based feature selection.
        keep &= np.nan_to_num(np.nanstd(x[train], axis=0), nan=0) > 0
        imp = SimpleImputer(strategy="median")
        clean = np.empty((len(x), int(keep.sum())))
        clean[train] = imp.fit_transform(x[train][:, keep])
        clean[~train] = imp.transform(x[~train][:, keep])
        joblib.dump(imp, out / f"{ep}-imputer.joblib")
        write_json(out / f"{ep}-columns.json", [c for c, k in zip(columns, keep, strict=True) if k])
        if cfg["variant"] == "stack_mordred":
            features = pl.read_parquet(next(inp.glob(f"*-{ep}-features.parquet")))
            joined = frame.select("compound_id", "fold").join(
                features, on="compound_id", how="left", validate="1:1"
            )
            assert (joined["fold"] == joined["fold_right"]).all()
            base = np.asarray(joined["features"].to_list(), dtype=np.float32)
            assert base.shape == (len(x), 260)
            clean = np.column_stack((base, clean))
        binary = ep.endswith("is_TDI")
        kind = "classifier" if binary else "regressor"
        cls = TabPFNClassifier if binary else TabPFNRegressor
        model = cls.create_default_for_version(
            ModelVersion.V2_6,
            model_path=str(next(inp.glob(f"*-tabpfn-v2.6-{kind}-v2.6_default.ckpt"))),
            device="cuda",
            n_estimators=8,
            random_state=42 + fold,
            ignore_pretraining_limits=True,
            fit_mode="fit_preprocessors",
            n_preprocessing_jobs=1,
            memory_saving_mode="auto",
            show_progress_bar=False,
        )
        write_json(out / f"{ep}-parameters.json", model.get_params(deep=False))
        model.fit(clean[train], y[train])
        p = model.predict_proba(clean[~train])[:, 1] if binary else model.predict(clean[~train])
        assert np.isfinite(p).all()
        rows = (
            frame.filter(pl.Series(~train))
            .select("compound_id")
            .with_columns(
                pl.lit(ep).alias("endpoint"),
                pl.lit(fold).alias("fold"),
                pl.lit(cfg["split_run"]).alias("split_id"),
                pl.lit("umap").alias("strategy"),
                pl.Series("y_true", y[~train]),
                pl.Series("prediction", p),
            )
        )
        all_rows.append(rows)
        score = dict(
            endpoint=ep,
            fold=fold,
            n=len(p),
            **metrics(
                y[~train],
                p,
                frame["conf_low"].to_numpy()[~train],
                frame["conf_high"].to_numpy()[~train],
                binary,
            ),
        )
        scores.append(score)
        model.save_fit_state(out / f"{ep}.tabpfn_fit")
        print(json.dumps(score), flush=True)
        del model
        gc.collect()
        torch.cuda.empty_cache()
    pl.concat(all_rows).write_parquet(out / "oof.parquet")
    write_json(out / "summary.json", dict(config=cfg, metrics=scores))


if __name__ == "__main__":
    main()
