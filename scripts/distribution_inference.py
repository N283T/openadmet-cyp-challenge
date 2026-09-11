"""Recover mean/median/mode/quantiles from unchanged regression fold models."""

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


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    scratch = Path(os.environ["TMPDIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    torch.set_num_threads(4)
    aliases = (
        pl.read_parquet(inp / cfg["aliases"])
        .filter(pl.col("role") == "test")
        .select("compound_id", pl.col("molecule_name").alias("Molecule_Name"))
        .unique()
    )
    ref = pl.read_csv(inp / cfg["reference"])
    test = (
        ref.select("SMILES", "Molecule_Name")
        .join(aliases, on="Molecule_Name", how="left", validate="1:1")
        .sort("compound_id")
    )
    assert test.height == 750 and test["compound_id"].null_count() == 0
    foundation = pl.concat([pl.read_parquet(inp / p) for p in cfg["foundation"]])
    desc = pl.concat([pl.read_parquet(inp / p) for p in cfg["mordred"]])
    base = test.select("compound_id").join(foundation, on="compound_id", how="left", validate="1:1")
    raw = (
        base.select(sorted(c for c in base.columns if c.startswith("chemeleon_")))
        .to_numpy()
        .astype(np.float32)
    )
    all_predictions = []
    all_oof = []
    audits = []
    for item in cfg["folds"]:
        fold = item["fold"]
        auxiliary = np.zeros((750, 4), dtype=np.float32)
        for parent in item["auxiliary"]:
            meta = json.loads((inp / parent["summary"]).read_text())
            assert fold in meta["config"]["excluded_folds"]
            pred = pl.concat([pl.read_parquet(inp / p) for p in parent["predictions"]])
            values = (
                test.select("compound_id")
                .join(pred, on="compound_id", how="left", validate="1:1")
                .select(meta["heads"])
                .to_numpy()
            )
            auxiliary += values / 4
        original = pl.read_parquet(inp / item["oof"])
        for ep, paths in item["endpoints"].items():
            binary = ep.endswith("is_TDI")
            kind = "classifier" if binary else "regressor"
            # Rewrite only the copied archive's base-weight path, keeping originals untouched.
            remapped = scratch / f"{fold}-{ep}.tabpfn_fit"
            with (
                zipfile.ZipFile(inp / paths["model"]) as source,
                zipfile.ZipFile(remapped, "w") as dest,
            ):
                for name in source.namelist():
                    payload = source.read(name)
                    if name == "init_params.json":
                        params = json.loads(payload)
                        params["model_path"] = str(inp / cfg["weights"][kind])
                        payload = json.dumps(params).encode()
                    dest.writestr(name, payload)
            model = load_fitted_tabpfn_model(remapped, device="cuda")
            columns = json.loads((inp / paths["columns"]).read_text())
            imputer = joblib.load(inp / paths["imputer"])
            pca = joblib.load(inp / paths["pca"])

            def descriptors(ids, columns=columns, imputer=imputer):
                values = (
                    ids.join(desc, on="compound_id", how="left", validate="1:1")
                    .select(columns)
                    .to_numpy()
                    .astype(np.float64)
                )
                values[~np.isfinite(values)] = np.nan
                return imputer.transform(values)

            held = (
                pl.read_parquet(inp / paths["validation_features"])
                .filter(pl.col("fold") == fold)
                .sort("compound_id")
            )
            held_x = np.column_stack(
                (
                    np.asarray(held["features"].to_list(), dtype=np.float32),
                    descriptors(held.select("compound_id")),
                )
            )
            assert not binary
            held_dist = model.predict(held_x, output_type="main", quantiles=cfg["quantiles"])
            reloaded = held_dist["mean"]
            expected = (
                held.select("compound_id")
                .join(original.filter(pl.col("endpoint") == ep), on="compound_id", validate="1:1")[
                    "prediction"
                ]
                .to_numpy()
            )
            np.testing.assert_allclose(reloaded, expected, rtol=1e-5, atol=1e-5)
            audits.append(
                dict(
                    fold=fold,
                    endpoint=ep,
                    max_abs_difference=float(np.max(np.abs(reloaded - expected))),
                )
            )
            x = np.column_stack(
                (pca.transform(raw), auxiliary, descriptors(test.select("compound_id")))
            )
            dist = model.predict(x, output_type="main", quantiles=cfg["quantiles"])
            for population, values, frame, collector in [
                ("umap", held_dist, held.select("compound_id"), all_oof),
                (
                    "test",
                    dist,
                    test.select("compound_id", "Molecule_Name", "SMILES"),
                    all_predictions,
                ),
            ]:
                points = {k: values[k] for k in ["mean", "median", "mode"]}
                points.update(
                    {
                        f"q{q:.2f}": v
                        for q, v in zip(cfg["quantiles"], values["quantiles"], strict=True)
                    }
                )
                for estimator, pred in points.items():
                    assert np.isfinite(pred).all()
                    collector.append(
                        frame.with_columns(
                            pl.lit(fold).alias("fold"),
                            pl.lit(ep).alias("endpoint"),
                            pl.lit(population).alias("split_id"),
                            pl.lit(estimator).alias("estimator"),
                            pl.Series("prediction", np.asarray(pred, dtype=float)),
                        )
                    )
            print(json.dumps(audits[-1]), flush=True)
            del model
            gc.collect()
            torch.cuda.empty_cache()
    predictions = pl.concat(all_predictions)
    predictions.write_parquet(out / "fold_predictions.parquet")
    pl.concat(all_oof).write_parquet(out / "distribution-oof.parquet")
    write_json(
        out / "summary.json",
        dict(
            config=cfg,
            reload_audit=audits,
            scope="Unchanged TabPFN 2.6 regression model; multiple distribution summaries",
        ),
    )


if __name__ == "__main__":
    main()
