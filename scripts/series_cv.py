"""Shared structure-series CV and matched descriptor learner comparison."""

import gc
import importlib.util
import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
import torch
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.ML.Cluster import Butina
from scipy.stats import spearmanr
from sklearn.impute import SimpleImputer
from tabpfn import TabPFNRegressor

from cyp_challenge.splits import grouped_folds
from cyp_challenge.storage import write_json

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(name):
    return pl.read_parquet(next(inp.glob("*-" + name + ".parquet")))


comp, aliases, labels = read("compounds"), read("aliases"), read("measurements")
test = set(aliases.filter(pl.col("role") == "test")["compound_id"])
comp = comp.filter(~pl.col("compound_id").is_in(list(test))).sort("compound_id")
ids = comp["compound_id"].to_list()
mols = [Chem.MolFromSmiles(s) for s in comp["standardized_smiles"]]
# Collapse stereoisomers to connectivity blocks before clustering.
keys = [Chem.MolToInchiKey(m).split("-")[0] for m in mols]
unique = sorted(set(keys))
representative = {k: keys.index(k) for k in unique}
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fps = [gen.GetFingerprint(mols[representative[k]]) for k in unique]
dist = []
for i in range(1, len(fps)):
    dist.extend(1 - np.asarray(DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])))
clusters = Butina.ClusterData(
    dist, len(fps), 1 - cfg["similarity_threshold"], isDistData=True, reordering=True
)
keygroup = {unique[i]: g for g, cluster in enumerate(clusters) for i in cluster}
groups = [str(keygroup[k]) for k in keys]
folds = grouped_folds(groups, 5)
assign = comp.select("compound_id").with_columns(
    pl.Series("group", groups), pl.Series("fold", folds), pl.lit("series").alias("strategy")
)
seeds = []
for enzyme in ["CYP1A2", "CYP2C9", "CYP3A4"]:
    f = (
        labels.filter(pl.col("endpoint") == enzyme + "_pIC50_direct_inhibition")
        .join(assign, on="compound_id", validate="1:1")
        .sort(["value", "compound_id"], descending=[True, False])
    )
    chosen = f.unique(subset="group", keep="first", maintain_order=True).head(
        cfg["seeds_per_enzyme"]
    )
    for r in chosen.iter_rows(named=True):
        seeds.append(
            dict(enzyme=enzyme, compound_id=r["compound_id"], group=r["group"], pic50=r["value"])
        )
seedtable = pl.DataFrame(seeds)
seedtable.write_csv(out / "hit-seeds.csv")
origin = {
    g: ",".join(sorted(set(seedtable.filter(pl.col("group") == g)["enzyme"])))
    for g in seedtable["group"].unique()
}
assign = assign.with_columns(pl.Series("hit_origins", [origin.get(g, "") for g in groups]))
assign.write_parquet(out / "assignments.parquet")
old = (
    read("assignments")
    .filter(pl.col("strategy") == "umap")
    .select("compound_id", "fold", "strategy")
    .join(assign.select("compound_id", "hit_origins"), on="compound_id", validate="1:1")
)
allassign = pl.concat([assign.select(old.columns), old])
reg = labels.filter(pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition"))
coverage = reg.join(allassign, on="compound_id", validate="m:m").with_columns(
    (pl.col("hit_origins") != "").alias("hit_series")
)
coverage.group_by("strategy", "fold", "endpoint", "hit_series").agg(
    pl.len().alias("n"),
    pl.col("value").mean().alias("mean"),
    pl.col("value").median().alias("median"),
    (pl.col("value") < 4).sum().alias("below4"),
).write_csv(out / "coverage.csv")
desc = pl.concat([pl.read_parquet(p) for p in sorted(inp.glob("*-mordred-*.parquet"))])
cols = [c for c in desc.columns if c != "compound_id"]
spec = importlib.util.spec_from_file_location(
    "official", next(inp.glob("*-custom_scoring_functions.py"))
)
official = importlib.util.module_from_spec(spec)
spec.loader.exec_module(official)
torch.set_num_threads(4)
rows = []
if list(inp.glob("*-oof.parquet")):
    previous = read("oof").filter(pl.col("strategy") == "umap").drop("hit_origins")
    previous = previous.join(
        assign.select("compound_id", "hit_origins"), on="compound_id", validate="m:1"
    )
    order = [
        "compound_id",
        "endpoint",
        "fold",
        "strategy",
        "hit_origins",
        "y_true",
        "conf_low",
        "conf_high",
        "model",
        "prediction",
    ]
    rows.append(previous.select(order))
for strategy in cfg.get("strategies", ["umap", "series"]):
    for ep in sorted(reg["endpoint"].unique()):
        f = (
            reg.filter(pl.col("endpoint") == ep)
            .join(
                allassign.filter(pl.col("strategy") == strategy), on="compound_id", validate="1:1"
            )
            .join(desc, on="compound_id", validate="1:1")
            .sort("compound_id")
        )
        x = f.select(cols).to_numpy().astype(float)
        x[~np.isfinite(x)] = np.nan
        y = f["value"].to_numpy()
        for fold in range(5):
            va = f["fold"].to_numpy() == fold
            tr = ~va
            imputer = SimpleImputer(strategy="median")
            xx = imputer.fit_transform(x[tr])
            vx = imputer.transform(x[va])
            for name in ["lightgbm", "tabpfn"]:
                if name == "lightgbm":
                    model = lgb.LGBMRegressor(
                        objective="regression_l1",
                        n_estimators=300,
                        num_leaves=15,
                        min_child_samples=20,
                        learning_rate=0.05,
                        reg_lambda=1,
                        n_jobs=4,
                        random_state=42 + fold,
                        verbosity=-1,
                        deterministic=True,
                        force_col_wise=True,
                    )
                else:
                    model = TabPFNRegressor(
                        n_estimators=8,
                        device="cuda",
                        random_state=42 + fold,
                        model_path=str(next(inp.glob("*.ckpt"))),
                        ignore_pretraining_limits=True,
                    )
                model.fit(xx, y[tr])
                p = model.predict(vx)
                assert np.isfinite(p).all()
                rows.append(
                    f.filter(pl.Series(va))
                    .select(
                        "compound_id",
                        "endpoint",
                        "fold",
                        "strategy",
                        "hit_origins",
                        pl.col("value").alias("y_true"),
                        "conf_low",
                        "conf_high",
                    )
                    .with_columns(
                        pl.lit(name).alias("model"), pl.Series("prediction", p.astype(float))
                    )
                )
                print(
                    json.dumps(
                        dict(strategy=strategy, endpoint=ep, fold=fold, model=name, n=int(va.sum()))
                    ),
                    flush=True,
                )
                del model
                gc.collect()
                torch.cuda.empty_cache()
        pl.concat(rows).write_parquet(out / f"progress-{strategy}-{ep}.parquet")
oof = pl.concat(rows)
oof.write_parquet(out / "oof.parquet")
metrics = []
for (strategy, ep, model), f in oof.group_by("strategy", "endpoint", "model"):
    for subset in ["all", "hit_series"]:
        g = f if subset == "all" else f.filter(pl.col("hit_origins") != "")
        for fold in [-1, 0, 1, 2, 3, 4]:
            h = g if fold == -1 else g.filter(pl.col("fold") == fold)
            if h.height < 2:
                continue
            y, p, lo, hi = [
                h[c].to_numpy() for c in ["y_true", "prediction", "conf_low", "conf_high"]
            ]
            metrics.append(
                dict(
                    strategy=strategy,
                    endpoint=ep,
                    model=model,
                    subset=subset,
                    fold=fold,
                    n=h.height,
                    ST_RAE=float(
                        official.rae_soft_threshold_absolute_error(
                            y, p, y_true_lower=lo, y_true_upper=hi
                        )
                    ),
                    MAE=float(np.abs(y - p).mean()),
                    Spearman=float(spearmanr(y, p).statistic),
                )
            )
m = pl.DataFrame(metrics)
m.write_csv(out / "metrics.csv")
summary = (
    m.filter(pl.col("fold") == -1)
    .group_by("strategy", "model", "subset")
    .agg(pl.col("ST_RAE", "MAE", "Spearman").mean())
)
summary.write_csv(out / "summary.csv")
write_json(
    out / "summary.json",
    dict(
        config=cfg,
        groups=len(clusters),
        compounds=len(ids),
        hit_groups=len(origin),
        results=summary.to_dicts(),
        scope=(
            "Structure-only folds; outcome-ranked seeds define descriptive subsets. "
            "Seeds held out with groups. Training-only imputation; no auxiliary labels."
        ),
    ),
)
print(summary)
