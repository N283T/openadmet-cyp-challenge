"""Matched challenge CV using only the frozen scratch-PubChem600 embedding."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
from tabpfn import TabPFNRegressor
from tabpfn.constants import ModelVersion

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(n):
    return pl.read_parquet(next(inp.glob("*-" + n + ".parquet")))


a, labels, emb = read("parent-assignments"), read("measurements"), read("pubchem600-embeddings")
cols = sorted(c for c in emb.columns if c.startswith("emb_"))
assert len(cols) == 600
torch.set_num_threads(4)
rows = []
for ep in sorted(
    labels.filter(pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition"))["endpoint"].unique()
):
    h = (
        labels.filter(pl.col("endpoint") == ep)
        .join(a, on="compound_id", validate="1:1")
        .join(emb, on="compound_id", validate="1:1")
        .sort("compound_id")
    )
    x, y = h.select(cols).to_numpy(), h["value"].to_numpy()
    for fold in range(5):
        va = (h["fold"].to_numpy() == fold) & ~h["is_anchor"].to_numpy()
        model = TabPFNRegressor.create_default_for_version(
            ModelVersion.V2_6,
            n_estimators=8,
            model_path=str(inp / cfg["weights"]),
            device="cuda",
            random_state=42 + fold,
            auto_scale_n_estimators=False,
            n_preprocessing_jobs=1,
            memory_saving_mode="auto",
            show_progress_bar=False,
            ignore_pretraining_limits=True,
            fit_mode="fit_preprocessors",
        )
        model.fit(x[~va], y[~va])
        pred = model.predict(x[va]).astype(float)
        assert np.isfinite(pred).all()
        model.save_fit_state(out / f"{ep}-fold{fold}.tabpfn_fit")
        rows.append(
            h.filter(pl.Series(va))
            .select(
                "compound_id",
                "endpoint",
                "fold",
                "is_neighbor",
                "conf_low",
                "conf_high",
                pl.col("value").alias("y_true"),
            )
            .with_columns(
                pl.lit("pubchem600_tabpfn").alias("method"),
                pl.lit(cfg["split_run"]).alias("split_id"),
                pl.Series("prediction", pred),
            )
        )
        print(ep, fold, int(va.sum()), flush=True)
pl.concat(rows).write_parquet(out / "pubchem600-oof.parquet")
(out / "summary.json").write_text(
    json.dumps(dict(config=cfg, models=20, features=600, feature_additions=False), indent=2)
)
