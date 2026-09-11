"""Inspect pretrained scores as features across all train and blinded test molecules."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import pearsonr, rankdata, spearmanr

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])


def read(name):
    return pl.read_parquet(next(inp.glob("*-" + name + ".parquet")))


raw = read("direct-predictions")
base = raw.select("compound_id", "fold").unique().sort("compound_id")
for score in ["probability", "logit"]:
    wide = raw.with_columns(
        (pl.col("endpoint").str.split("_").list.first() + "_" + score).alias("feature")
    ).pivot(on="feature", index="compound_id", values=score)
    base = base.join(wide, on="compound_id", validate="1:1", maintain_order="left")
scores = sorted(c for c in base.columns if c.endswith(("_probability", "_logit")))
heads = sorted(c for c in read("oof-log2fc").columns if "@" in c)
aux = pl.concat(
    [
        read("oof-log2fc").select("compound_id", *heads),
        read("test-log2fc").select("compound_id", *heads),
    ]
)
base = base.join(aux, on="compound_id", validate="1:1", maintain_order="left").with_columns(
    pl.when(pl.col("fold") == -1)
    .then(pl.lit("test"))
    .otherwise(pl.lit("train_all"))
    .alias("population")
)
base.write_parquet(out / "score-features.parquet")
assert base.height == 6896 and base["compound_id"].n_unique() == 6896
summary, rows, distributions = [], [], []
for pop in ["train_all", "test"]:
    h = base.filter(pl.col("population") == pop)
    for a in scores:
        x = h[a].to_numpy()
        distributions.append(
            dict(
                population=pop,
                feature=a,
                n=len(x),
                mean=float(x.mean()),
                sd=float(x.std()),
                q05=float(np.quantile(x, 0.05)),
                median=float(np.median(x)),
                q95=float(np.quantile(x, 0.95)),
            )
        )
        for b in heads + [s for s in scores if s.endswith("_probability")]:
            rows.append(
                dict(
                    population=pop,
                    a=a,
                    b=b,
                    n=h.height,
                    pearson=float(pearsonr(x, h[b]).statistic),
                    spearman=float(spearmanr(x, h[b]).statistic),
                )
            )
pl.DataFrame(rows).write_csv(out / "score-correlations.csv")
pl.DataFrame(distributions).write_csv(out / "score-distributions.csv")
labels = read("measurements").filter(pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition"))
label_rows = []
for (ep,), labels_ep in labels.group_by("endpoint"):
    h = labels_ep.join(base, on="compound_id", validate="1:1")
    for s in scores:
        label_rows.append(
            dict(
                endpoint=ep,
                feature=s,
                n=h.height,
                pearson=float(pearsonr(h[s], h["value"]).statistic),
                spearman=float(spearmanr(h[s], h["value"]).statistic),
            )
        )
pl.DataFrame(label_rows).write_csv(out / "all-label-correlations.csv")


def correlations(x, y):
    xc = x - x.mean(axis=0)
    yc = y - y.mean(axis=0)
    return (xc.T @ yc) / (np.linalg.norm(xc, axis=0)[:, None] * np.linalg.norm(yc, axis=0)[None, :])


feature_rows = []
for family in ["mordred", "chemeleon"]:
    feature = pl.concat(pl.read_parquet(p) for p in sorted(inp.glob("*-" + family + "-*.parquet")))
    names = [c for c in feature.columns if c != "compound_id"]
    h = base.join(feature, on="compound_id", validate="1:1", maintain_order="left")
    for pop in ["train_all", "test"]:
        g = h.filter(pl.col("population") == pop)
        x = g.select(scores).to_numpy().astype(float)
        y = g.select(names).to_numpy().astype(float)
        complete = np.isfinite(y).all(axis=0) & (np.nanstd(y, axis=0) > 0)
        js = np.flatnonzero(complete)
        pc = correlations(x, y[:, js])
        sc = correlations(rankdata(x, axis=0), rankdata(y[:, js], axis=0))
        for j, colidx in enumerate(js):
            for i, s in enumerate(scores):
                feature_rows.append(
                    dict(
                        population=pop,
                        family=family,
                        score=s,
                        feature=names[colidx],
                        n=g.height,
                        pearson=float(pc[i, j]),
                        spearman=float(sc[i, j]),
                    )
                )
        for j in np.flatnonzero(~complete):
            valid = np.isfinite(y[:, j])
            if valid.sum() < 30 or np.std(y[valid, j]) == 0:
                continue
            pc = correlations(x[valid], y[valid, j, None])
            sc = correlations(rankdata(x[valid], axis=0), rankdata(y[valid, j, None], axis=0))
            for i, s in enumerate(scores):
                feature_rows.append(
                    dict(
                        population=pop,
                        family=family,
                        score=s,
                        feature=names[j],
                        n=int(valid.sum()),
                        pearson=float(pc[i, 0]),
                        spearman=float(sc[i, 0]),
                    )
                )
        print(f"{family} {pop} complete", flush=True)
features = pl.DataFrame(feature_rows)
features.write_parquet(out / "all-feature-correlations.parquet")
top = (
    features.with_columns(pl.col("spearman").abs().alias("abs_spearman"))
    .sort("abs_spearman", descending=True)
    .group_by("population", "family", "score", maintain_order=True)
    .head(5)
)
top.write_csv(out / "top-feature-correlations.csv")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            train_compounds=6146,
            test_compounds=750,
            admetica_features=scores,
            pred_log2fc_train="four eligible models excluding recipient fold",
            pred_log2fc_test="ten-model mean",
            feature_correlation_rows=features.height,
            missing_descriptors="pairwise finite observations, at least 30; constants omitted",
            scope="Feature correlations use no labels; test targets unavailable.",
        ),
        indent=2,
    )
)
print(
    pl.DataFrame(rows).filter(
        pl.col("a").str.ends_with("_probability") & pl.col("b").str.contains("@")
    )
)
