"""Infer blinded test predictions from immutable fold models and verify OOF reload."""

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
from cyp_challenge.tracks import ENDPOINTS


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
            reloaded = model.predict_proba(held_x)[:, 1] if binary else model.predict(held_x)
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
            pred = model.predict_proba(x)[:, 1] if binary else model.predict(x)
            assert np.isfinite(pred).all()
            all_predictions.append(
                test.select("compound_id", "Molecule_Name", "SMILES").with_columns(
                    pl.lit(fold).alias("fold"),
                    pl.lit(ep).alias("endpoint"),
                    pl.Series("prediction", pred),
                )
            )
            print(json.dumps(audits[-1]), flush=True)
            del model
            gc.collect()
            torch.cuda.empty_cache()
    predictions = pl.concat(all_predictions)
    predictions.write_parquet(out / "fold_predictions.parquet")
    mean = predictions.group_by("compound_id", "Molecule_Name", "SMILES", "endpoint").agg(
        pl.col("prediction").mean()
    )
    mean.write_parquet(out / "ensemble_probabilities.parquet")
    wide = mean.pivot(
        on="endpoint", index=["compound_id", "Molecule_Name", "SMILES"], values="prediction"
    ).sort("Molecule_Name")
    wide.select("SMILES", "Molecule_Name", *ENDPOINTS["regression"]).write_parquet(
        out / "regression.parquet"
    )
    wide.select(
        "SMILES",
        "Molecule_Name",
        *[(pl.col(ep) >= cfg["threshold"]).alias(ep) for ep in ENDPOINTS["classification"]],
    ).write_parquet(out / "classification.parquet")
    write_json(
        out / "summary.json",
        dict(
            config=cfg,
            reload_audit=audits,
            classification_positive_counts={
                ep: int((wide[ep] >= cfg["threshold"]).sum()) for ep in ENDPOINTS["classification"]
            },
        ),
    )


if __name__ == "__main__":
    main()
