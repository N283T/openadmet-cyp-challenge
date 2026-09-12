"""Pair-excluded primary-screen LGBM features, with membership and routing audits."""

import itertools
import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from scipy.stats import spearmanr

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(i):
    return pl.read_parquet(next(inp.glob(f"{i:03d}-*")))


a, targets, emb, scores, old = [read(i) for i in range(5)]
print("schemas", a.schema, targets.schema, scores.columns[:10], flush=True)
sc = [c for c in scores.columns if c != "compound_id" and scores[c].dtype.is_numeric()]
ec = [c for c in emb.columns if c.startswith("emb_")]
features = (
    emb.select("compound_id", *ec)
    .join(scores.select("compound_id", *sc), on="compound_id", validate="1:1")
    .sort("compound_id")
)
assert features.height == 6896 and len(ec) == 600 and len(sc) == 47
cols = ec + sc
x = features.select(cols).to_numpy()
assert np.isfinite(x).all()
h = features.select("compound_id").join(
    a, on="compound_id", how="left", validate="1:1", maintain_order="left"
)
heads = sorted(targets["head"].unique())
assert len(heads) == 4
assert not targets.select("compound_id", "head").is_duplicated().any()
t = h.join(
    targets.pivot(on="head", index="compound_id", values="value"),
    on="compound_id",
    how="left",
    validate="1:1",
    maintain_order="left",
)
y = t.select(heads).to_numpy()
fold = t["fold"].fill_null(-1).to_numpy()
anchor = t["is_anchor"].fill_null(False).to_numpy()
assert a.height == 6146
assert a.group_by("group").agg(pl.col("fold").n_unique().alias("n"))["n"].max() == 1
models = {}
audits = []
for pair in itertools.combinations(range(5), 2):
    excluded = np.isin(fold, pair) & ~anchor
    pred = np.zeros((len(x), 4))
    for j, head in enumerate(heads):
        train = (fold >= 0) & ~excluded & np.isfinite(y[:, j])
        assert not (train & excluded).any()
        model = lgb.LGBMRegressor(**cfg["lgbm"])
        model.fit(x[train], y[train, j])
        pred[:, j] = model.predict(x)
        model.booster_.save_model(str(out / f"pair{pair[0]}{pair[1]}-head{j}.txt"))
        t.select("compound_id", "fold", "is_anchor").with_columns(
            pl.Series("has_label", np.isfinite(y[:, j])),
            pl.Series("is_training", train),
            pl.Series("is_excluded", excluded),
        ).write_parquet(out / f"pair{pair[0]}{pair[1]}-head{j}-membership.parquet")
        audits.append(
            dict(
                pair=list(pair),
                head=head,
                training=int(train.sum()),
                excluded_labels=int((excluded & np.isfinite(y[:, j])).sum()),
                overlap=0,
            )
        )
    assert np.isfinite(pred).all()
    models[pair] = pred
    features.select("compound_id").with_columns(
        *[pl.Series(head, pred[:, j]) for j, head in enumerate(heads)]
    ).write_parquet(out / f"pair{pair[0]}{pair[1]}-predictions.parquet")
    print("completed", pair, flush=True)
frames = []
for outer in range(5):
    associated = [p for p in models if outer in p]
    vals = np.stack([models[p] for p in associated]).mean(axis=0)
    for recipient in range(5):
        if recipient == outer:
            continue
        mask = (fold == recipient) & ~anchor
        vals[mask] = models[tuple(sorted((outer, recipient)))][mask]
    for index in np.flatnonzero(fold >= 0):
        used = (
            associated
            if (anchor[index] or fold[index] == outer)
            else [tuple(sorted((outer, int(fold[index]))))]
        )
        assert np.allclose(
            vals[index], np.stack([models[p][index] for p in used]).mean(axis=0), atol=1e-12
        )
        assert all(outer in p for p in used)
        if not anchor[index]:
            assert all(int(fold[index]) in p for p in used)
    frames.append(
        t.select("compound_id", "fold", "is_anchor").with_columns(
            pl.lit(outer).alias("outer_fold"),
            pl.Series("aux_model_count", np.where((fold == outer) | anchor | (fold < 0), 4, 1)),
            *[pl.Series(f"aux_{j}", vals[:, j]) for j in range(4)],
        )
    )
allnested = pl.concat(frames)
allnested.write_parquet(out / "all-nested-log2fc.parquet")
nested = allnested.filter(pl.col("fold").is_not_null())
assert nested.height == 30730
nested.write_parquet(out / "nested-log2fc.parquet")
rows = []
for j, head in enumerate(heads):
    measured = (
        targets.filter(pl.col("head") == head)
        .join(a.select("compound_id", "fold", "is_anchor"), on="compound_id")
        .filter(~pl.col("is_anchor"))
    )
    for name, source in [("new_lgbm", nested), ("old_chemprop", old)]:
        z = measured.join(
            source.filter(pl.col("outer_fold") == pl.col("fold")).select("compound_id", f"aux_{j}"),
            on="compound_id",
            validate="1:1",
        )
        yy, pp = z["value"].to_numpy(), z[f"aux_{j}"].to_numpy()
        rows.append(
            dict(
                head=head,
                method=name,
                n=z.height,
                rmse=float(np.sqrt(np.mean((yy - pp) ** 2))),
                mae=float(np.mean(np.abs(yy - pp))),
                pearson=float(np.corrcoef(yy, pp)[0, 1]),
                spearman=float(spearmanr(yy, pp).statistic),
            )
        )
pl.DataFrame(rows).write_csv(out / "primary-metrics.csv")
(out / "membership-audit.json").write_text(json.dumps(audits, indent=2))
(out / "summary.json").write_text(
    json.dumps(
        dict(
            heads=heads,
            features=len(cols),
            rows=nested.height,
            all_rows=allnested.height,
            anchor_policy=(
                "Anchors retained in primary fitting and average four outer-excluding models. "
                "All nonanchor outer and recipient folds excluded."
            ),
            config=cfg,
        ),
        indent=2,
    )
)
print(pl.DataFrame(rows), flush=True)
