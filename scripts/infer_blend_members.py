"""Reload compact and auxiliary-pretrained TabPFN members for OOF/test distributions."""

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
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    fold = cfg["fold"]
    torch.set_num_threads(4)

    def table(name):
        return pl.read_parquet(inp / name)

    test = (
        table(cfg["aliases"])
        .filter(pl.col("role") == "test")
        .select("compound_id")
        .unique()
        .sort("compound_id")
    )
    assert test.height == 750
    oofrows, testrows, audit = [], [], []
    for variant in ["compact", "pretrained"]:
        item = cfg[variant]
        old = table(item["oof"])
        if variant == "compact":
            foundation = pl.concat([table(p) for p in cfg["foundation"]])
            raw = (
                test.join(foundation, on="compound_id", validate="1:1")
                .select(sorted(c for c in foundation.columns if c.startswith("chemeleon_")))
                .to_numpy()
                .astype(np.float32)
            )
            auxiliary = np.zeros((750, 4), dtype=np.float32)
            for source in cfg["auxiliary"]:
                meta = json.loads((inp / source["summary"]).read_text())
                assert fold in meta["config"]["excluded_folds"]
                data = pl.concat([table(p) for p in source["predictions"]])
                auxiliary += (
                    test.join(data, on="compound_id", validate="1:1")
                    .select(meta["heads"])
                    .to_numpy()
                    / 4
                )
        else:
            emb = pl.concat([table(p) for p in item["embeddings"]])
            assert set(emb["fold"]) == {fold}
            emb = emb.drop("fold")
            testx = np.asarray(
                test.join(emb, on="compound_id", validate="1:1")["embedding"].to_list(),
                dtype=np.float32,
            )
        for ep, paths in item["endpoints"].items():
            held = old.filter(pl.col("endpoint") == ep).sort("compound_id")
            if variant == "compact":
                features = table(paths["features"]).filter(pl.col("fold") == fold)
                heldx = np.asarray(
                    held.select("compound_id")
                    .join(features, on="compound_id", validate="1:1")["features"]
                    .to_list(),
                    dtype=np.float32,
                )
                pca = joblib.load(inp / paths["pca"])
                testx = np.column_stack([pca.transform(raw), auxiliary])
            else:
                heldx = np.asarray(
                    held.select("compound_id")
                    .join(emb, on="compound_id", validate="1:1")["embedding"]
                    .to_list(),
                    dtype=np.float32,
                )
            remapped = Path(os.environ["TMPDIR"]) / f"{variant}-{ep}.tabpfn_fit"
            with (
                zipfile.ZipFile(inp / paths["model"]) as source,
                zipfile.ZipFile(remapped, "w") as dest,
            ):
                for name in source.namelist():
                    payload = source.read(name)
                    if name == "init_params.json":
                        params = json.loads(payload)
                        params["model_path"] = str(inp / cfg["weights"])
                        payload = json.dumps(params).encode()
                    dest.writestr(name, payload)
            model = load_fitted_tabpfn_model(remapped, device="cuda")
            for split, x, frame, collector in [
                ("umap", heldx, held.select("compound_id"), oofrows),
                ("test", testx, test, testrows),
            ]:
                values = model.predict(x, output_type="main", quantiles=cfg["quantiles"])
                if split == "umap":
                    np.testing.assert_allclose(
                        values["mean"], held["prediction"], rtol=1e-5, atol=1e-5
                    )
                    audit.append(
                        dict(
                            variant=variant,
                            endpoint=ep,
                            max_abs_diff=float(
                                np.max(np.abs(values["mean"] - held["prediction"].to_numpy()))
                            ),
                        )
                    )
                points = {k: values[k] for k in ["mean", "median"]}
                points.update(
                    {
                        f"q{q:.2f}": v
                        for q, v in zip(cfg["quantiles"], values["quantiles"], strict=True)
                    }
                )
                for estimator, p in points.items():
                    assert np.isfinite(p).all()
                    collector.append(
                        frame.with_columns(
                            pl.lit(ep).alias("endpoint"),
                            pl.lit(fold).alias("fold"),
                            pl.lit(split).alias("split_id"),
                            pl.lit(variant).alias("member"),
                            pl.lit(estimator).alias("estimator"),
                            pl.Series("prediction", np.asarray(p, dtype=float)),
                        )
                    )
            print(json.dumps(audit[-1]), flush=True)
            del model
            gc.collect()
            torch.cuda.empty_cache()
    pl.concat(oofrows).write_parquet(out / "member-oof.parquet")
    pl.concat(testrows).write_parquet(out / "member-test.parquet")
    write_json(out / "summary.json", dict(config=cfg, audit_passed=True, reload_audit=audit))


if __name__ == "__main__":
    main()
