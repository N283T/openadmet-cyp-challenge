"""CYP2D6 AID1851 efficacy-only external trees with explicit challenge exclusion."""

import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupShuffleSplit

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(i):
    return pl.read_parquet(next(inp.glob(f"{i:03d}-*")))


raw, challenge = read(0), read(1).sort("compound_id")
raw = raw.filter(pl.col("isoform") == "CYP2D6")
forbidden = {
    Chem.MolToInchiKey(Chem.MolFromSmiles(s)).split("-")[0]
    for s in challenge["standardized_smiles"]
}
assert (
    not raw.filter(~pl.col("challenge_overlap") & pl.col("connectivity").is_not_null())[
        "connectivity"
    ]
    .is_in(list(forbidden))
    .any()
)
accepted = raw.filter(~pl.col("challenge_overlap") & pl.col("connectivity").is_not_null())
continuous = ["Max_Response", "Activity at 57.14 uM", "Activity at 1.022 uM"]
accepted = accepted.with_columns(
    *[pl.col(c).cast(pl.Float64, strict=False) for c in continuous],
    pl.when(pl.col("PUBCHEM_ACTIVITY_OUTCOME") == "Inactive")
    .then(1.0)
    .when(pl.col("PUBCHEM_ACTIVITY_OUTCOME") == "Active")
    .then(0.0)
    .otherwise(None)
    .alias("inactive"),
)
summary = accepted.group_by("PUBCHEM_ACTIVITY_OUTCOME").agg(
    pl.len().alias("rows"), *[pl.col(c).is_finite().sum().alias(c) for c in continuous]
)
summary.write_csv(out / "coverage.csv")
curated = (
    accepted.group_by("compound_id", "standardized_smiles", "connectivity")
    .agg(
        *[pl.col(c).filter(pl.col(c).is_finite()).median() for c in continuous],
        pl.when(pl.col("inactive").drop_nulls().n_unique() == 1)
        .then(pl.col("inactive").drop_nulls().first())
        .otherwise(None)
        .alias("inactive"),
    )
    .sort("compound_id")
)
curated.write_parquet(out / "efficacy-targets.parquet")
fp = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def matrix(ss):
    return np.asarray([fp.GetFingerprintAsNumPy(Chem.MolFromSmiles(s)) for s in ss])


x = matrix(curated["standardized_smiles"])
xc = matrix(challenge["standardized_smiles"])
groups = [
    MurckoScaffold.MurckoScaffoldSmiles(smiles=s) or c
    for s, c in curated.select("standardized_smiles", "connectivity").iter_rows()
]
tr, va = next(GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=42).split(x, groups=groups))
isval = np.isin(np.arange(len(x)), va)
curated.select("compound_id", "connectivity").with_columns(
    pl.Series("validation", isval), pl.Series("scaffold", groups)
).write_parquet(out / "external-memberships.parquet")
features = challenge.select("compound_id")
metrics = []
holdouts = []
for name, column, objective, clip in [
    ("max_l1", "Max_Response", "regression_l1", False),
    ("max_l2", "Max_Response", "regression", True),
    ("high_dose", "Activity at 57.14 uM", "regression_l1", False),
    ("mid_dose", "Activity at 1.022 uM", "regression_l1", False),
    ("inactive", "inactive", "binary", False),
]:
    yy = curated[column].to_numpy().astype(float)
    valid = np.isfinite(yy)
    if clip:
        yy = np.clip(yy, -150, 50)
    train = valid & ~isval
    val = valid & isval
    params = dict(
        objective=objective,
        n_estimators=1000,
        num_leaves=31,
        min_child_samples=40,
        learning_rate=0.04,
        reg_lambda=5,
        n_jobs=3,
        random_state=42,
        verbosity=-1,
        deterministic=True,
        force_col_wise=True,
    )
    cls = lgb.LGBMClassifier if objective == "binary" else lgb.LGBMRegressor
    model = cls(**params).fit(
        x[train],
        yy[train],
        eval_set=[(x[val], yy[val])],
        callbacks=[lgb.early_stopping(75, verbose=False)],
    )
    pv = model.predict_proba(x[val])[:, 1] if objective == "binary" else model.predict(x[val])
    stats = dict(
        method=name,
        n=int(valid.sum()),
        validation=int(val.sum()),
        trees=model.best_iteration_,
        spearman=float(spearmanr(yy[val], pv).statistic),
        MAE=float(np.mean(np.abs(yy[val] - pv))),
        constant_MAE=float(np.mean(np.abs(yy[val] - np.median(yy[train])))),
    )
    if objective == "binary":
        stats["ROC_AUC"] = float(roc_auc_score(yy[val], pv))
    metrics.append(stats)
    print(stats, flush=True)
    holdouts.append(
        curated.filter(pl.Series(val))
        .select("compound_id")
        .with_columns(
            pl.lit(name).alias("method"), pl.Series("y", yy[val]), pl.Series("prediction", pv)
        )
    )
    model.booster_.save_model(str(out / f"{name}-holdout.lgbm.txt"))
    params["n_estimators"] = model.best_iteration_
    final = cls(**params).fit(x[valid], yy[valid])
    final.booster_.save_model(str(out / f"{name}-full.lgbm.txt"))
    pc = final.predict_proba(xc)[:, 1] if objective == "binary" else final.predict(xc)
    features = features.with_columns(pl.Series("efficacy_" + name, pc))
features.write_parquet(out / "efficacy-features.parquet")
pl.concat(holdouts).write_parquet(out / "external-holdout.parquet")
pl.DataFrame(metrics).write_csv(out / "external-metrics.csv")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            config=cfg,
            raw_rows=raw.height,
            accepted_rows=accepted.height,
            unique_compounds=curated.height,
            challenge_connectivity_excluded=True,
            metrics=metrics,
            meaning="External efficacy and inactivity probabilities; not pIC50",
        ),
        indent=2,
    )
)
