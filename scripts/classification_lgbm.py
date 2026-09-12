"""Fixed UMAP TDI screen with inner-only MCC thresholds and frozen external features."""

import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(name):
    return pl.read_parquet(inp / cfg[name])


comp = read("compounds").sort("compound_id")
ids = comp.select("compound_id")
labels = read("labels").filter(pl.col("endpoint").is_in(["CYP2D6_is_TDI", "CYP3A4_is_TDI"]))
assert len(labels) == 5081
folds = (
    read("assignments").filter(pl.col("strategy") == "umap").select("compound_id", "fold", "group")
)
test = read("test_mapping").select("compound_id").unique().sort("compound_id")
assert len(test) == 750
index = {c: i for i, c in enumerate(ids["compound_id"])}
itest = np.array([index[c] for c in test["compound_id"]])


def matrix(table, prefix=None):
    cols = [
        c for c in table.columns if c != "compound_id" and (prefix is None or c.startswith(prefix))
    ]
    result = (
        ids.join(
            table.select("compound_id", *cols),
            on="compound_id",
            how="left",
            validate="1:1",
            maintain_order="left",
        )
        .select(cols)
        .to_numpy()
        .astype(np.float32)
    )
    result[~np.isfinite(result)] = np.nan
    return result


mordred = pl.concat([pl.read_parquet(inp / name) for name in cfg["mordred"]])
xm = matrix(mordred)
xc = matrix(read("class_emb"), "emb_")
xmix = matrix(read("mixed_emb"), "emb_")
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
xf = np.stack(
    [gen.GetFingerprintAsNumPy(Chem.MolFromSmiles(s)) for s in comp["standardized_smiles"]]
).astype(np.float32)
views = {
    "morgan": xf,
    "mordred": xm,
    "fp_mordred": np.column_stack([xf, xm]),
    "external_class": xc,
    "external_mixed": xmix,
    "external_class_mordred": np.column_stack([xc, xm]),
}


def threshold(y, p):
    order = np.argsort(-p, kind="stable")
    yy = y[order]
    pp = p[order]
    k = np.r_[np.flatnonzero(pp[:-1] != pp[1:]) + 1, len(p)]
    tp = np.cumsum(yy)[k - 1]
    fp = k - tp
    fn = y.sum() - tp
    tn = len(y) - y.sum() - fp
    den = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    scores = np.divide(tp * tn - fp * fn, den, out=np.zeros_like(den, dtype=float), where=den > 0)
    thresholds = pp[k - 1]
    choices = [
        (float(s), -abs(float(t) - 0.5), float(t)) for s, t in zip(scores, thresholds, strict=True)
    ]
    choices.append((float(matthews_corrcoef(y, p >= 0.5)), 0.0, 0.5))
    return max(choices)[2]


def metric(y, p, binary):
    return dict(
        MCC=float(matthews_corrcoef(y, binary)),
        AUROC=float(roc_auc_score(y, p)),
        AP=float(average_precision_score(y, p)),
        positives=int(binary.sum()),
        n=len(y),
    )


rows, tests, stats, configs, memberships = [], [], [], [], []
for (ep,), lab in labels.group_by("endpoint"):
    h = lab.join(folds, on="compound_id", validate="1:1").sort("compound_id")
    assert len(h) == len(lab)
    yy = h["value"].to_numpy().astype(int)
    ff = h["fold"].to_numpy()
    ix = np.array([index[c] for c in h["compound_id"]])
    for name, x in views.items():
        for balanced in [False, True]:
            method = name + ("_balanced" if balanced else "_plain")
            for fold in range(5):
                train = ff != fold
                val = ~train
                inner = ff == (fold + 1) % 5
                base = dict(
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
                    class_weight="balanced" if balanced else None,
                )
                provisional = lgb.LGBMClassifier(**base).fit(
                    x[ix[train & ~inner]], yy[train & ~inner]
                )
                innerp = provisional.predict_proba(x[ix[inner]])[:, 1]
                cut = threshold(yy[inner], innerp)
                provisional.booster_.save_model(str(out / f"{method}-{ep}-{fold}-inner.txt"))
                model = lgb.LGBMClassifier(**base).fit(x[ix[train]], yy[train])
                p = model.predict_proba(x[ix[val]])[:, 1]
                tp = model.predict_proba(x[itest])[:, 1]
                assert np.isfinite(p).all() and np.isfinite(tp).all()
                model.booster_.save_model(str(out / f"{method}-{ep}-{fold}.txt"))
                pred = (
                    h.filter(pl.Series(val))
                    .select(
                        "compound_id", "endpoint", "fold", "group", pl.col("value").alias("y_true")
                    )
                    .with_columns(
                        pl.Series("prediction", p),
                        pl.lit(method).alias("method"),
                        pl.lit(cut).alias("threshold"),
                        pl.lit(cfg["split_run"]).alias("split_id"),
                    )
                )
                rows.append(pred)
                tests.append(
                    test.with_columns(
                        pl.lit(ep).alias("endpoint"),
                        pl.lit(fold).alias("fold"),
                        pl.lit(method).alias("method"),
                        pl.Series("prediction", tp),
                        pl.lit(cut).alias("threshold"),
                        pl.lit(cfg["split_run"]).alias("split_id"),
                    )
                )
                for policy, tt in [("fixed_05", 0.5), ("inner_mcc", cut)]:
                    stats.append(
                        dict(
                            endpoint=ep,
                            method=method,
                            fold=fold,
                            policy=policy,
                            **metric(yy[val], p, p >= tt),
                        )
                    )
                configs.append(
                    dict(
                        endpoint=ep,
                        method=method,
                        fold=fold,
                        threshold=cut,
                        n_train=int(train.sum()),
                        n_inner_fit=int((train & ~inner).sum()),
                        n_inner_validation=int(inner.sum()),
                        dimensions=x.shape[1],
                        params=base,
                    )
                )
                if name == "morgan" and not balanced:
                    memberships.append(
                        h.select("compound_id", "endpoint", "fold").with_columns(
                            pl.lit(fold).alias("outer_fold"),
                            pl.Series("outer_validation", val),
                            pl.Series("inner_validation", inner),
                        )
                    )
                print(ep, method, fold, flush=True)
result = pl.concat(rows)
result.write_parquet(out / "oof.parquet")
pl.concat(tests).write_parquet(out / "fold-test-probabilities.parquet")
for (ep, method), h in result.group_by("endpoint", "method"):
    y = h["y_true"].to_numpy().astype(int)
    p = h["prediction"].to_numpy()
    for policy, cut in [("fixed_05", 0.5), ("inner_mcc", h["threshold"].to_numpy())]:
        stats.append(
            dict(endpoint=ep, method=method, fold=-1, policy=policy, **metric(y, p, p >= cut))
        )
pl.DataFrame(stats).write_csv(out / "metrics.csv")
pl.concat(memberships).write_parquet(out / "membership.parquet")
(out / "models.json").write_text(json.dumps(configs, indent=2))
(out / "summary.json").write_text(
    json.dumps(
        dict(
            methods=12,
            folds=5,
            models=len(configs),
            features="Frozen public encoders and static molecular descriptors only",
            threshold="Fixed inner UMAP fold; outer fold excluded from fit and threshold choice",
            split_run=cfg["split_run"],
        ),
        indent=2,
    )
)
print(
    pl.DataFrame(stats)
    .filter(pl.col("fold") == -1)
    .sort("endpoint", "MCC", descending=[False, True])
    .write_csv(),
    flush=True,
)
