"""Matched-feature TabICLv2 point regression with fixed distribution readouts."""

import gc
import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
import torch

from cyp_challenge.storage import write_json


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
        from tabicl import TabICLRegressor

        model = TabICLRegressor(
            n_estimators=8,
            batch_size=1,
            device="cuda",
            random_state=42 + fold,
            n_jobs=4,
            model_path=str(inp / cfg["tabicl_weights"]),
            allow_auto_download=False,
            use_amp="auto",
            use_fa3="auto",
            offload_mode="auto",
            disk_offload_dir=os.environ["TMPDIR"],
        )
        model.fit(x[tr], frame["value"].to_numpy()[tr])
        for split, xx, ids, collector in [
            ("umap", x[va], frame.filter(pl.Series(va)).select("compound_id"), rows),
            ("test", tx, test, ts),
        ]:
            dist = model.predict(xx, output_type=["mean", "median", "quantiles"], alphas=[0.55])
            for method, pred in [
                ("mean", dist["mean"]),
                ("median", dist["median"]),
                ("q55", np.asarray(dist["quantiles"]).reshape(-1)),
            ]:
                assert len(pred) == ids.height and np.isfinite(pred).all()
                collector.append(
                    ids.with_columns(
                        pl.lit(ep).alias("endpoint"),
                        pl.lit(fold).alias("fold"),
                        pl.lit(split).alias("split_id"),
                        pl.lit("tabicl_" + method).alias("method"),
                        pl.Series("prediction", np.asarray(pred, dtype=float)),
                    )
                )
        info.append(
            dict(
                endpoint=ep,
                fold=fold,
                n_train=int(tr.sum()),
                n_features=x.shape[1],
                parameters={k: str(v) for k, v in model.get_params().items()},
            )
        )
        print(json.dumps(info[-1]), flush=True)
        del model
        gc.collect()
        torch.cuda.empty_cache()
    pl.concat(rows).write_parquet(out / "tabicl-oof.parquet")
    pl.concat(ts).write_parquet(out / "tabicl-test.parquet")
    write_json(
        out / "summary.json",
        dict(
            config=cfg,
            audit_passed=True,
            targets=info,
            scope="Matched features and UMAP folds; wrapper defaults, fixed TabICLv2 weights",
        ),
    )


if __name__ == "__main__":
    main()
