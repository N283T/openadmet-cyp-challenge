"""Controlled fixed-validation TDI learning size and shared-population comparison."""

import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.metrics import average_precision_score, roc_auc_score

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(name):
    return pl.read_parquet(inp / cfg[name])


comp = read("compounds").sort("compound_id")
assign = (
    read("assignments").filter(pl.col("strategy") == "umap").select("compound_id", "fold", "group")
)
labels = read("labels").filter(pl.col("endpoint").is_in(["CYP2D6_is_TDI", "CYP3A4_is_TDI"]))
mordred = pl.concat([pl.read_parquet(inp / name) for name in cfg["mordred"]])
cols = [c for c in mordred.columns if c != "compound_id"]
xm = (
    comp.select("compound_id")
    .join(mordred, on="compound_id", how="left", validate="1:1", maintain_order="left")
    .select(cols)
    .to_numpy()
    .astype(np.float32)
)
xm[~np.isfinite(xm)] = np.nan
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
xf = np.stack(
    [gen.GetFingerprintAsNumPy(Chem.MolFromSmiles(s)) for s in comp["standardized_smiles"]]
).astype(np.float32)
x = np.column_stack([xf, xm])
index = {c: i for i, c in enumerate(comp["compound_id"])}
frames = {
    ep: labels.filter(pl.col("endpoint") == ep)
    .join(assign, on="compound_id", validate="1:1")
    .sort("compound_id")
    for ep in ["CYP2D6_is_TDI", "CYP3A4_is_TDI"]
}
shared = set(frames["CYP2D6_is_TDI"]["compound_id"]) & set(frames["CYP3A4_is_TDI"]["compound_id"])
reference = read("reference").filter(pl.col("method") == "fp_mordred_balanced")
predictions, memberships, metrics = [], [], []


def stratified(frame, n, seed):
    assert 0 < n <= frame.height
    groups = frame.group_by("fold", "value", maintain_order=True).agg(pl.col("compound_id"))
    sizes = np.array([len(ids) for ids in groups["compound_id"]])
    quota = n * sizes / sizes.sum()
    allocation = np.floor(quota).astype(int)
    remainder = n - allocation.sum()
    allocation[np.argsort(-(quota - allocation), kind="stable")[:remainder]] += 1
    rng = np.random.default_rng(seed)
    chosen = []
    for ids, count in zip(groups["compound_id"], allocation, strict=True):
        chosen.extend(rng.choice(ids.to_numpy(), size=int(count), replace=False).tolist())
    assert len(chosen) == n and len(set(chosen)) == n
    return frame.filter(pl.col("compound_id").is_in(chosen))


def fit_eval(method, ep, fold, train, val, existing=None):
    ti = np.array([index[c] for c in train["compound_id"]])
    vi = np.array([index[c] for c in val["compound_id"]])
    y = train["value"].to_numpy().astype(int)
    if existing is None:
        model = lgb.LGBMClassifier(
            objective="binary",
            n_estimators=400,
            num_leaves=15,
            min_child_samples=20,
            learning_rate=0.03,
            reg_lambda=2,
            n_jobs=3,
            random_state=42 + fold,
            verbosity=-1,
            deterministic=True,
            force_col_wise=True,
            class_weight="balanced",
        ).fit(x[ti], y)
        booster = model.booster_
        booster.save_model(str(out / f"{method}-{ep}-{fold}.txt"))
    else:
        booster = lgb.Booster(model_file=str(inp / existing))
    for scope, data, ii in [("train", train, ti), ("validation", val, vi)]:
        p = booster.predict(x[ii], num_threads=3)
        yy = data["value"].to_numpy().astype(int)
        metrics.append(
            dict(
                method=method,
                endpoint=ep,
                fold=fold,
                scope=scope,
                n=len(yy),
                positive_rate=float(yy.mean()),
                AUROC=roc_auc_score(yy, p),
                AP=average_precision_score(yy, p),
            )
        )
        predictions.append(
            data.select(
                "compound_id", "endpoint", "fold", "group", pl.col("value").alias("y_true")
            ).with_columns(
                pl.lit(fold).alias("outer_fold"),
                pl.lit(method).alias("method"),
                pl.lit(scope).alias("scope"),
                pl.Series("prediction", p),
            )
        )
        if existing and scope == "validation":
            expected = reference.filter((pl.col("endpoint") == ep) & (pl.col("fold") == fold))
            pair = (
                data.select("compound_id")
                .with_columns(pl.Series("p", p))
                .join(
                    expected.select("compound_id", "prediction"), on="compound_id", validate="1:1"
                )
            )
            assert pair.height == data.height
            np.testing.assert_allclose(pair["p"], pair["prediction"], rtol=0, atol=1e-12)
    memberships.append(
        train.select(
            "compound_id", "endpoint", "fold", pl.col("value").alias("y_true")
        ).with_columns(pl.lit(method).alias("method"), pl.lit(fold).alias("outer_fold"))
    )
    pl.DataFrame(metrics).write_csv(out / "fold-metrics.csv")
    pl.concat(predictions).write_parquet(out / "predictions.parquet")
    pl.concat(memberships).write_parquet(out / "selected-training-ids.parquet")
    print(method, ep, fold, train.height, val.height, flush=True)


for fold in range(5):
    if "route_ids" in cfg:
        allowed = read("route_ids").filter(
            (pl.col(cfg["route_column"]) == cfg["route_value"])
            & (pl.col("endpoint") == "CYP3A4_is_TDI")
        )["compound_id"]
        a4 = frames["CYP3A4_is_TDI"].filter(pl.col("compound_id").is_in(allowed.implode()))
        assert a4.height == 2334
        train, val = a4.filter(pl.col("fold") != fold), a4.filter(pl.col("fold") == fold)
        fit_eval("paired_assay_full", "CYP3A4_is_TDI", fold, train, val)
        n = frames["CYP2D6_is_TDI"].filter(pl.col("fold") != fold).height
        for seed in cfg["seeds"]:
            fit_eval(
                f"paired_assay_d6_size_seed{seed}",
                "CYP3A4_is_TDI",
                fold,
                stratified(train, n, seed),
                val,
            )
        continue
    for ep, frame in frames.items():
        fit_eval(
            "full",
            ep,
            fold,
            frame.filter(pl.col("fold") != fold),
            frame.filter(pl.col("fold") == fold),
            cfg["full_models"][f"{ep}-{fold}"],
        )
    a4 = frames["CYP3A4_is_TDI"]
    train, val = a4.filter(pl.col("fold") != fold), a4.filter(pl.col("fold") == fold)
    n = frames["CYP2D6_is_TDI"].filter(pl.col("fold") != fold).height
    for seed in cfg["seeds"]:
        matched = stratified(train, n, seed)
        half = stratified(matched, n // 2, seed + 1000)
        fit_eval(f"d6_size_seed{seed}", "CYP3A4_is_TDI", fold, matched, val)
        fit_eval(f"half_d6_size_seed{seed}", "CYP3A4_is_TDI", fold, half, val)
    for ep, frame in frames.items():
        paired = frame.filter(pl.col("compound_id").is_in(list(shared)))
        fit_eval(
            "shared_population",
            ep,
            fold,
            paired.filter(pl.col("fold") != fold),
            paired.filter(pl.col("fold") == fold),
        )
validation = pl.concat(predictions).filter(pl.col("scope") == "validation")
summary = []
for (method, ep), table in validation.group_by("method", "endpoint"):
    summary.append(
        dict(
            method=method,
            endpoint=ep,
            n=table.height,
            AUROC=roc_auc_score(table["y_true"], table["prediction"]),
            AP=average_precision_score(table["y_true"], table["prediction"]),
        )
    )
pl.DataFrame(summary).write_csv(out / "pooled-metrics.csv")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            new_fits=20 if "route_ids" in cfg else 40,
            reused_full_models=0 if "route_ids" in cfg else 10,
            shared_compounds=len(shared),
            seeds=cfg["seeds"],
            original_full_predictions_exact=None if "route_ids" in cfg else True,
        ),
        indent=2,
    )
)
