"""Replace nested primary features with fixed previously selected FT readout parameters."""

import importlib.util
import itertools
import json
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def path(name):
    return inp / cfg["inputs"][name]


def read(name):
    return pl.read_parquet(path(name))


a, labels, aux = read("assignments"), read("labels"), read("aux")
labels = labels.filter(pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition"))
cols = [f"emb_{i:03d}" for i in range(600)] + [f"aux_{i}" for i in range(4)]
checks, rows, fits = [], [], []
# Verify all primary pair model training memberships and exact nested routing.
pair_predictions = {}
for pair in itertools.combinations(range(5), 2):
    name = "".join(map(str, pair))
    p = read("pair" + name).join(a.select("compound_id"), on="compound_id").sort("compound_id")
    pc = [c for c in p.columns if c != "compound_id"]
    assert len(pc) == 4
    pair_predictions[pair] = p.select(pc).to_numpy()
    for head in range(4):
        m = read(f"pair{name}-head{head}")
        expected = m["fold"].is_in(pair).fill_null(False) & ~m["is_anchor"].fill_null(False)
        assert m["is_excluded"].equals(expected.alias("is_excluded"))
        assert m.filter(pl.col("is_training") & pl.col("is_excluded")).is_empty()
for fold in range(5):
    heldout = a.filter((pl.col("fold") == fold) & ~pl.col("is_anchor"))
    m = read(f"membership{fold}")
    assert set(m.filter(pl.col("role") == "outer_holdout")["compound_id"]) == set(
        heldout["compound_id"]
    )
    assert (
        m.filter(pl.col("role").is_in(["training", "inner_validation"]))
        .join(heldout, on="compound_id")
        .is_empty()
    )
    na = aux.filter(pl.col("outer_fold") == fold).sort("compound_id")
    aa = a.sort("compound_id")
    assert na["compound_id"].equals(aa["compound_id"])
    pp = [pair for pair in pair_predictions if fold in pair]
    expected = np.stack([pair_predictions[pair] for pair in pp]).mean(0)
    for recipient in range(5):
        if recipient == fold:
            continue
        mask = ((aa["fold"] == recipient) & ~aa["is_anchor"]).to_numpy()
        expected[mask] = pair_predictions[tuple(sorted((fold, recipient)))][mask]
    assert np.allclose(
        expected, na.select([f"aux_{i}" for i in range(4)]).to_numpy(), atol=1e-12, rtol=0
    )
    emb = read(f"emb{fold}")
    assert emb["fold"].unique().to_list() == [fold]
    features = a.join(
        emb.select("compound_id", *cols[:600]), on="compound_id", validate="1:1"
    ).join(na.select("compound_id", *cols[600:]), on="compound_id", validate="1:1")
    for ep in sorted(labels["endpoint"].unique()):
        h = (
            labels.filter(pl.col("endpoint") == ep)
            .join(features, on="compound_id", validate="1:1")
            .sort("compound_id")
        )
        va = ((h["fold"] == fold) & ~h["is_anchor"]).to_numpy()
        x, y = h.select(cols).to_numpy().astype(np.float32), h["value"].to_numpy()
        assert np.isfinite(x).all()
        selected = json.loads(path(f"params{fold}-{ep}").read_text())
        params = dict(
            selected["params"],
            n_jobs=3,
            random_state=42 + fold,
            verbosity=-1,
            deterministic=True,
            force_col_wise=True,
        )
        model = lgb.LGBMRegressor(**params)
        model.fit(x[~va], y[~va])
        prediction = model.predict(x[va])
        assert np.isfinite(prediction).all()
        model.booster_.save_model(str(out / f"{ep}-{fold}.txt"))
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
                pl.Series("prediction", prediction),
                pl.lit(cfg["method"]).alias("method"),
                pl.lit(cfg["split_run"]).alias("split_id"),
            )
        )
        fits.append(
            dict(
                fold=fold,
                endpoint=ep,
                train=int((~va).sum()),
                validation=int(va.sum()),
                params=params,
                original_trial=selected["trial"],
            )
        )
        print(json.dumps(fits[-1]), flush=True)
    checks.append(dict(fold=fold, outer_exclusion=True, primary_routing_exact=True))
oof = pl.concat(rows)
keys = ["compound_id", "endpoint", "fold"]
assert oof.height == 6384 and oof.select(keys).unique().height == 6384
oof.write_parquet(out / "external-oof.parquet")
old, final = read("old").sort(keys), read("final").sort(keys)
oof = oof.sort(keys)
for f in [old, final]:
    assert f.select(keys + ["y_true", "conf_low", "conf_high"]).equals(
        oof.select(keys + ["y_true", "conf_low", "conf_high"])
    )
aligned = oof.with_columns(
    old["prediction"].alias("old"),
    final["prediction"].alias("final"),
    pl.col("prediction").alias("new"),
).with_columns(
    (0.75 * pl.col("final") + 0.25 * pl.col("new")).alias("add25"),
    ((pl.col("final") + pl.col("new")) / 2).alias("add50"),
    (pl.col("final") + (pl.col("new") - pl.col("old")) / 5).alias("replace_component"),
)
spec = importlib.util.spec_from_file_location("official", path("scorer"))
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
metrics = []
for subset in ["all", "neighbors"]:
    frame = aligned if subset == "all" else aligned.filter(pl.col("is_neighbor"))
    for (ep,), g in frame.group_by("endpoint"):
        for method in ["old", "new", "final", "add25", "add50", "replace_component"]:
            score = scorer.rae_soft_threshold_absolute_error(
                g["y_true"].to_numpy(),
                g[method].to_numpy(),
                y_true_lower=g["conf_low"].to_numpy(),
                y_true_upper=g["conf_high"].to_numpy(),
            )
            metrics.append(
                dict(endpoint=ep, subset=subset, method=method, n=g.height, ST_RAE=score)
            )
metrics = pl.DataFrame(metrics).sort("subset", "method", "endpoint")
metrics.write_csv(out / "comparison.csv")
metrics.group_by("subset", "method").agg(pl.col("ST_RAE").mean()).sort(
    "subset", "ST_RAE"
).write_csv(out / "macro.csv")
aligned.write_parquet(out / "aligned-oof.parquet")
(out / "summary.json").write_text(
    json.dumps(dict(fits=fits, checks=checks, rows=6384, retuned=False, config=cfg), indent=2)
)
