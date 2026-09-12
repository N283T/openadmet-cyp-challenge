"""Paired fixed-CV substrate additions; control exactly matches the primary readout."""

import json
import os
from pathlib import Path

import joblib
import lightgbm as lgb
import numpy as np
import polars as pl
import torch
from tabpfn import TabPFNRegressor
from tabpfn.constants import ModelVersion

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
torch.set_num_threads(4)


def read(index):
    return pl.read_parquet(next(inp.glob(f"{index:03d}-*")))


emb, labels, assignments, nested, substrate, previous = [read(i) for i in range(6)]
ecols = sorted(c for c in emb.columns if c.startswith("emb_"))
assert emb["compound_id"].n_unique() == substrate["compound_id"].n_unique() == 6896
assert set(emb["compound_id"]) == set(substrate["compound_id"])
rows, fits, memberships = [], [], []
for ep in sorted(
    labels.filter(pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition"))["endpoint"].unique()
):
    h = (
        labels.filter(pl.col("endpoint") == ep)
        .join(assignments, on="compound_id", validate="1:1")
        .join(emb, on="compound_id", validate="1:1")
        .join(substrate, on="compound_id", validate="1:1")
        .sort("compound_id")
    )
    y = h["value"].to_numpy()
    for fold in range(5):
        va = (h["fold"].to_numpy() == fold) & ~h["is_anchor"].to_numpy()
        auxiliary = h.select("compound_id").join(
            nested.filter(pl.col("outer_fold") == fold),
            on="compound_id",
            how="left",
            validate="1:1",
            maintain_order="left",
        )
        x = np.column_stack(
            [
                h.select(ecols).to_numpy(),
                auxiliary.select([f"aux_{i}" for i in range(4)]).to_numpy(),
            ]
        )
        assert np.isfinite(x).all()
        memberships.append(
            h.select("compound_id", "is_anchor").with_columns(
                pl.lit(ep).alias("endpoint"),
                pl.lit(fold).alias("outer_fold"),
                pl.Series("validation", va),
            )
        )
        for method, jcols in cfg["arms"].items():
            xx = x
            if jcols:
                raw = h.select(jcols).to_numpy().astype(float)
                finite = np.isfinite(raw)
                keep = np.array(
                    [np.unique(raw[~va & finite[:, j], j]).size > 1 for j in range(len(jcols))]
                )
                kept = raw[:, keep]
                missing = ~np.isfinite(kept)
                medians = np.array(
                    [np.median(kept[~va & ~missing[:, j], j]) for j in range(kept.shape[1])]
                )
                imputed = np.where(missing, medians, kept)
                # Always append an indicator for retained descriptors: a validation-only
                # failure remains representable without inspecting validation statistics.
                xx = np.column_stack([x, imputed, missing.astype(float)])
                joblib.dump(
                    dict(columns=jcols, keep=keep, medians=medians),
                    out / f"{method}-{ep}-{fold}-preprocessor.joblib",
                )
            params = dict(
                objective="regression_l1",
                n_estimators=300,
                num_leaves=15,
                min_child_samples=20,
                learning_rate=0.05,
                reg_lambda=1,
                n_jobs=3,
                random_state=42 + fold,
                verbosity=-1,
                deterministic=True,
                force_col_wise=True,
            )
            if cfg.get("backend", "lgbm") == "lgbm":
                model = lgb.LGBMRegressor(**params)
            else:
                assert cfg["backend"] == "tabpfn"
                model = TabPFNRegressor.create_default_for_version(
                    ModelVersion.V2_6,
                    n_estimators=8,
                    model_path=str(next(inp.glob("006-*"))),
                    device="cuda",
                    random_state=42 + fold,
                    auto_scale_n_estimators=False,
                    n_preprocessing_jobs=1,
                    memory_saving_mode="auto",
                    show_progress_bar=False,
                    ignore_pretraining_limits=True,
                    fit_mode="fit_preprocessors",
                )
            model.fit(xx[~va], y[~va])
            pred = model.predict(xx[va]).astype(float)
            assert np.isfinite(pred).all()
            if cfg.get("backend", "lgbm") == "lgbm":
                model.booster_.save_model(str(out / f"{method}-{ep}-{fold}.lgbm.txt"))
            else:
                model.save_fit_state(out / f"{method}-{ep}-{fold}.tabpfn_fit")
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
                    pl.lit(method).alias("method"),
                    pl.lit(cfg["split_run"]).alias("split_id"),
                    pl.Series("prediction", pred),
                )
            )
            fits.append(
                dict(
                    method=method,
                    endpoint=ep,
                    fold=fold,
                    train=int((~va).sum()),
                    validation=int(va.sum()),
                    dimensions=xx.shape[1],
                )
            )
            print(json.dumps(fits[-1]), flush=True)
result = pl.concat(rows)
keys = ["compound_id", "endpoint", "fold"]
control = result.filter(pl.col("method") == "control").sort(keys)
previous = previous.sort(keys)
assert control.height == previous.height == 6384
for col in keys + ["y_true", "conf_low", "conf_high"]:
    assert control[col].equals(previous[col]), ("Control reproduction failed", col)
max_difference = float(
    np.max(np.abs(control["prediction"].to_numpy() - previous["prediction"].to_numpy()))
)
assert max_difference <= cfg.get("control_tolerance", 0), max_difference
result.write_parquet(out / "substrate-oof.parquet")
pl.concat(memberships).write_parquet(out / "memberships.parquet")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            config=cfg,
            embedding_columns=ecols,
            auxiliary_columns=[f"aux_{i}" for i in range(4)],
            effective_lgbm_parameters=params if cfg.get("backend", "lgbm") == "lgbm" else None,
            fits=fits,
            baseline_bitwise_equal=control["prediction"].equals(previous["prediction"]),
            baseline_max_abs_difference=max_difference,
            preprocessing="substrate mask and medians fit on outer training rows only",
        ),
        indent=2,
    )
)
