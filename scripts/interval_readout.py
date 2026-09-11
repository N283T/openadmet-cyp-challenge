"""Fit TabPFN to interval marginals; choose the median of their equal mixture.

For expected distance outside [L,U], derivative is F_L(q)+F_U(q)-1.
The median of an equal mixture of marginal L/U distributions is the Bayes action.
Quantile-grid approximation here uses no validation/test bounds as model inputs.
"""

import gc
import json
import os
import zipfile
from pathlib import Path

import joblib
import numpy as np
import polars as pl
import torch
from tabpfn.model_loading import load_fitted_tabpfn_model

from cyp_challenge.storage import write_json


def mixture_median(lower_samples, upper_samples):
    return np.median(np.column_stack([lower_samples, upper_samples]), axis=1)


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    torch.set_num_threads(4)

    def table(p):
        return pl.read_parquet(inp / p)

    test = (
        table(cfg["aliases"])
        .filter(pl.col("role") == "test")
        .select("compound_id")
        .unique()
        .sort("compound_id")
    )
    labels = table(cfg["measurements"])
    foundation = pl.concat([table(p) for p in cfg["foundation"]])
    desc = pl.concat([table(p) for p in cfg["mordred"]])
    raw = (
        test.join(foundation, on="compound_id", validate="1:1", maintain_order="left")
        .select(sorted(c for c in foundation.columns if c.startswith("chemeleon_")))
        .to_numpy()
        .astype(np.float32)
    )
    item = cfg["folds"][0]
    fold = item["fold"]
    aux = np.zeros((750, 4), dtype=np.float32)
    for parent in item["auxiliary"]:
        meta = json.loads((inp / parent["summary"]).read_text())
        assert fold in meta["config"]["excluded_folds"]
        data = pl.concat([table(p) for p in parent["predictions"]])
        aux += (
            test.join(data, on="compound_id", validate="1:1", maintain_order="left")
            .select(meta["heads"])
            .to_numpy()
            / 4
        )
    rows, ts, info = [], [], []
    for ep, paths in item["endpoints"].items():
        frame = table(paths["validation_features"]).sort("compound_id")
        frame = frame.join(
            labels.filter(pl.col("endpoint") == ep).select(
                "compound_id", "value", "conf_low", "conf_high"
            ),
            on="compound_id",
            validate="1:1",
            maintain_order="left",
        )
        assert frame["value"].null_count() == 0
        tr = frame["fold"].to_numpy() != fold
        va = ~tr
        imp = joblib.load(inp / paths["imputer"])
        cols = json.loads((inp / paths["columns"]).read_text())

        def descriptor(ids, cols=cols, imp=imp):
            x = (
                ids.join(desc, on="compound_id", validate="1:1", maintain_order="left")
                .select(cols)
                .to_numpy()
                .astype(float)
            )
            x[~np.isfinite(x)] = np.nan
            return imp.transform(x)

        x = np.column_stack(
            [
                np.asarray(frame["features"].to_list(), dtype=np.float32),
                descriptor(frame.select("compound_id")),
            ]
        )
        pca = joblib.load(inp / paths["pca"])
        tx = np.column_stack([pca.transform(raw), aux, descriptor(test)])
        if cfg.get("backend") == "lightgbm":
            import lightgbm as lgb

            y = frame["value"].to_numpy()
            lo, hi = frame["conf_low"].to_numpy(), frame["conf_high"].to_numpy()
            for target, trees in [("point", 300), ("interval", 300), ("interval", 800)]:
                train_x = x[tr] if target == "point" else np.vstack([x[tr], x[tr]])
                train_y = y[tr] if target == "point" else np.r_[lo[tr], hi[tr]]
                weights = np.ones(len(train_y)) if target == "point" else np.full(len(train_y), 0.5)
                model = lgb.LGBMRegressor(
                    objective="regression_l1",
                    n_estimators=trees,
                    learning_rate=0.05,
                    num_leaves=15,
                    min_child_samples=20,
                    reg_lambda=1.0,
                    n_jobs=4,
                    random_state=42 + fold,
                    deterministic=True,
                    force_col_wise=True,
                    verbosity=-1,
                )
                model.fit(train_x, train_y, sample_weight=weights)
                name = f"lgbm_{target}_{trees}"
                model.booster_.save_model(str(out / f"{ep}-{name}.txt"))
                for split, xx, ids, collector in [
                    ("umap", x[va], frame.filter(pl.Series(va)).select("compound_id"), rows),
                    ("test", tx, test, ts),
                ]:
                    pred = model.predict(xx)
                    assert np.isfinite(pred).all()
                    collector.append(
                        ids.with_columns(
                            pl.lit(ep).alias("endpoint"),
                            pl.lit(fold).alias("fold"),
                            pl.lit(split).alias("split_id"),
                            pl.lit(name).alias("method"),
                            pl.Series("prediction", pred),
                        )
                    )
                info.append(dict(endpoint=ep, method=name, n_train=int(tr.sum()), trees=trees))
                print(json.dumps(info[-1]), flush=True)
            continue
        archive = Path(os.environ["TMPDIR"]) / (ep + ".tabpfn_fit")
        with zipfile.ZipFile(inp / paths["model"]) as source, zipfile.ZipFile(archive, "w") as dest:
            for name in source.namelist():
                payload = source.read(name)
                if name == "init_params.json":
                    params = json.loads(payload)
                    params["model_path"] = str(inp / cfg["weights"]["regressor"])
                    payload = json.dumps(params).encode()
                dest.writestr(name, payload)
        model = load_fitted_tabpfn_model(archive, device="cuda")
        reference = table(item["oof"]).filter(pl.col("endpoint") == ep).sort("compound_id")
        np.testing.assert_allclose(
            model.predict(x[va]), reference["prediction"], rtol=1e-5, atol=1e-5
        )
        held, testdist = {}, {}
        targets = {
            "lower": frame["conf_low"].to_numpy(),
            "upper": frame["conf_high"].to_numpy(),
            "midpoint": (frame["conf_low"].to_numpy() + frame["conf_high"].to_numpy()) / 2,
        }
        for name, y in targets.items():
            assert np.isfinite(y).all()
            model.fit(x[tr], y[tr])
            hd = model.predict(x[va], output_type="main", quantiles=cfg["quantiles"])
            td = model.predict(tx, output_type="main", quantiles=cfg["quantiles"])
            held[name] = hd
            testdist[name] = td
            model.save_fit_state(out / f"{ep}-{name}.tabpfn_fit")
            info.append(
                dict(
                    endpoint=ep,
                    target=name,
                    fold=fold,
                    n_train=int(tr.sum()),
                    target_mean=float(y[tr].mean()),
                    target_sd=float(y[tr].std()),
                )
            )
            print(json.dumps(info[-1]), flush=True)
        for split, dist, ids, collector in [
            ("umap", held, frame.filter(pl.Series(va)).select("compound_id"), rows),
            ("test", testdist, test, ts),
        ]:
            ls, us = [np.column_stack(dist[k]["quantiles"]) for k in ["lower", "upper"]]
            # All endpoints' marginal samples retain their own learned shape.
            points = {
                "interval_mixture_median": mixture_median(ls, us),
                "midpoint_mean": dist["midpoint"]["mean"],
                "midpoint_median": dist["midpoint"]["median"],
                "average_bound_medians": (dist["lower"]["median"] + dist["upper"]["median"]) / 2,
            }
            sampleframe = ids.with_columns(
                pl.lit(ep).alias("endpoint"),
                pl.lit(fold).alias("fold"),
                pl.lit(split).alias("split_id"),
                pl.Series("lower_quantiles", ls),
                pl.Series("upper_quantiles", us),
            )
            sampleframe.write_parquet(out / f"{ep}-{split}-marginals.parquet")
            for name, p in points.items():
                assert np.isfinite(p).all()
                collector.append(
                    ids.with_columns(
                        pl.lit(ep).alias("endpoint"),
                        pl.lit(fold).alias("fold"),
                        pl.lit(split).alias("split_id"),
                        pl.lit(name).alias("method"),
                        pl.Series("prediction", np.asarray(p, dtype=float)),
                    )
                )
        del model
        gc.collect()
        torch.cuda.empty_cache()
    pl.concat(rows).write_parquet(out / "interval-oof.parquet")
    pl.concat(ts).write_parquet(out / "interval-test.parquet")
    write_json(
        out / "summary.json",
        dict(
            config=cfg,
            audit_passed=True,
            targets=info,
            scope=(
                "LightGBM point/interval L1 models fit on outer training only"
                if cfg.get("backend") == "lightgbm"
                else "Marginal interval models fit on outer training only; "
                "quantile-grid Bayes-action approximation"
            ),
        ),
    )


if __name__ == "__main__":
    main()
