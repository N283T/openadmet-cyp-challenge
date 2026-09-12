"""Fixed strong regularization and shuffled-training-label sanity controls."""

import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(name):
    return pl.read_parquet(inp / name)


comp = read(cfg["compounds"]).sort("compound_id")
labels = read(cfg["labels"]).filter(pl.col("endpoint") == "CYP2D6_is_TDI")
assigned = read(cfg["assignments"]).filter(pl.col("strategy") == "umap")
h = labels.join(assigned, on="compound_id", validate="1:1").sort("compound_id")
assert h.height == 1497
mordred = pl.concat([read(p) for p in cfg["mordred"]])
cols = [c for c in mordred.columns if c != "compound_id"]
xm = (
    comp.select("compound_id")
    .join(mordred, on="compound_id", how="left", validate="1:1", maintain_order="left")
    .select(cols)
    .to_numpy()
    .astype(np.float32)
)
xm[~np.isfinite(xm)] = np.nan
fg = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
xf = np.stack(
    [fg.GetFingerprintAsNumPy(Chem.MolFromSmiles(s)) for s in comp["standardized_smiles"]]
).astype(np.float32)
lookup = {c: i for i, c in enumerate(comp["compound_id"])}
ix = [lookup[c] for c in h["compound_id"]]
x = np.column_stack([xf, xm])[ix]
y = h["value"].to_numpy().astype(int)
folds = h["fold"].to_numpy()


def score(y, p):
    return dict(
        AUC=float(roc_auc_score(y, p)),
        AP=float(average_precision_score(y, p)),
        Brier=float(brier_score_loss(y, p)),
        n=len(y),
        prevalence=float(y.mean()),
    )


rows = []
preds = []
members = []
for recipe in ["strong_regularization", "train_label_permutation"]:
    for fold in range(5):
        train = folds != fold
        val = ~train
        seed = 42 + fold
        used = y[train].copy()
        if recipe == "train_label_permutation":
            used = np.random.default_rng(seed).permutation(used)
        params = dict(
            objective="binary",
            n_estimators=100 if recipe == "strong_regularization" else 400,
            num_leaves=7 if recipe == "strong_regularization" else 15,
            min_child_samples=80 if recipe == "strong_regularization" else 20,
            reg_lambda=10 if recipe == "strong_regularization" else 2,
            learning_rate=0.03,
            class_weight="balanced",
            n_jobs=3,
            random_state=seed,
            verbosity=-1,
            deterministic=True,
            force_col_wise=True,
        )
        model = lgb.LGBMClassifier(**params).fit(x[train], used)
        model.booster_.save_model(str(out / f"{recipe}-fold{fold}.txt"))
        tp = model.predict_proba(x[train])[:, 1]
        vp = model.predict_proba(x[val])[:, 1]
        for role, truth, p in [
            ("train_fit_labels", used, tp),
            ("train_original_labels", y[train], tp),
            ("heldout", y[val], vp),
        ]:
            rows.append(dict(method=recipe, fold=fold, role=role, **score(truth, p)))
        membership = h.select("compound_id", "endpoint", "fold", "group").with_columns(
            pl.lit(fold).alias("outer_fold"),
            pl.lit(seed).alias("seed"),
            pl.lit(recipe).alias("method"),
            pl.Series("training", train),
            pl.Series("original_label", y),
            pl.Series("fit_label", np.where(train, 0, -1)),
        )
        fitall = np.full(len(y), np.nan)
        fitall[train] = used
        members.append(membership.with_columns(pl.Series("fit_label", fitall)))
        preds.append(
            h.filter(pl.Series(val))
            .select("compound_id", "endpoint", "fold", "group")
            .with_columns(
                pl.Series("y_true", y[val]),
                pl.Series("probability", vp),
                pl.lit(recipe).alias("method"),
            )
        )
        trainpred = (
            h.filter(pl.Series(train))
            .select("compound_id", "endpoint", "fold", "group")
            .with_columns(
                pl.Series("y_true", y[train]),
                pl.Series("fit_label", used),
                pl.Series("probability", tp),
                pl.lit(recipe).alias("method"),
                pl.lit(fold).alias("outer_fold"),
            )
        )
        trainpred.write_parquet(out / f"{recipe}-fold{fold}-train.parquet")
        print(recipe, fold, score(y[val], vp), flush=True)
result = pl.concat(preds)
result.write_parquet(out / "oof.parquet")
pl.concat(members).write_parquet(out / "memberships.parquet")
for (method,), t in result.group_by("method"):
    rows.append(
        dict(
            method=method,
            fold=-1,
            role="heldout",
            **score(t["y_true"].to_numpy(), t["probability"].to_numpy()),
        )
    )
baseline = read(cfg["baseline"])
for method in ["fp_mordred_plain", "fp_mordred_balanced"]:
    t = baseline.filter(
        (pl.col("endpoint") == "CYP2D6_is_TDI") & (pl.col("method") == method)
    ).sort("compound_id")
    assert t.height == 1497 and t["compound_id"].equals(h["compound_id"])
    np.testing.assert_array_equal(t["y_true"], y)
    rows.append(
        dict(
            method="original_" + method,
            fold=-1,
            role="heldout",
            **score(y, t["prediction"].to_numpy()),
        )
    )
pl.DataFrame(rows).write_csv(out / "metrics.csv")
(out / "verification.json").write_text(
    json.dumps(
        dict(
            recipes=cfg["recipes"],
            fits=10,
            labels=1497,
            permutation="Only outer-training labels shuffled; true heldout labels unchanged; "
            "fitted and original train-label scores both saved",
            purpose="Sanity controls, not deployment selection",
            features=x.shape[1],
        ),
        indent=2,
    )
)
