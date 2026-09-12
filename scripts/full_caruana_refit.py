"""Refit the retained 24-candidate Caruana recipe on every endpoint label."""

import json
import os
from collections import defaultdict
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(key):
    return pl.read_parquet(inp / cfg[key])


compounds = read("compounds").sort("compound_id")
ids = compounds.select("compound_id")
labels = read("measurements").filter(pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition"))
test = read("test_mapping")
testids = test.select("compound_id").unique().sort("compound_id")
assert test.height == testids.height == 750
assert labels.height == labels.select("compound_id", "endpoint").n_unique() == 6525
assert labels.join(testids, on="compound_id", how="inner").is_empty()
selection = json.loads((inp / cfg["selection"]).read_text())
assert len(selection["selections"]) == 20
feature_columns = {}


def matrix(key, cols=None):
    table = read(key)
    cols = cols or sorted(c for c in table.columns if c != "compound_id")
    x = (
        ids.join(
            table.select("compound_id", *cols),
            on="compound_id",
            validate="1:1",
            maintain_order="left",
        )
        .select(cols)
        .to_numpy()
        .astype(np.float32)
    )
    assert np.isfinite(x).all(), key
    feature_columns[key] = cols
    return x


xc, xm, xt = matrix("class_emb"), matrix("mixed_emb"), matrix("tree")
xn = matrix("aux", [f"aux_{i}" for i in range(4)])
assert xc.shape[1] == xm.shape[1] == 600 and xt.shape[1] == 20
fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
xf = np.stack(
    [fpgen.GetFingerprintAsNumPy(Chem.MolFromSmiles(s)) for s in compounds["standardized_smiles"]]
)
views = dict(
    morgan=xf,
    logfc=xn,
    **{"class": xc},
    mixed=xm,
    tree=xt,
    class_logfc=np.column_stack([xc, xn]),
    mixed_logfc=np.column_stack([xm, xn]),
    tree_logfc=np.column_stack([xt, xn]),
)
index = {c: i for i, c in enumerate(ids["compound_id"])}
itest = np.array([index[c] for c in testids["compound_id"]])
rows, memberships, manifest = [], [], []
for ep in sorted(labels["endpoint"].unique()):
    selections = [s for s in selection["selections"] if s["endpoint"] == ep]
    assert {s["fold"] for s in selections} == set(range(5))
    weights, candidates = defaultdict(float), {}
    for selected in selections:
        assert np.isclose(sum(c["weight"] for c in selected["selected"]), 1)
        for c in selected["selected"]:
            weights[c["name"]] += c["weight"] / 5
            params = dict(c["params"], random_state=42)
            if c["name"] in candidates:
                assert candidates[c["name"]]["params"] == params
            candidates[c["name"]] = {"view": c["view"], "params": params}
    assert np.isclose(sum(weights.values()), 1)
    train = labels.filter(pl.col("endpoint") == ep).sort("compound_id")
    itr = np.array([index[c] for c in train["compound_id"]])
    y = train["value"].to_numpy()
    assert np.isfinite(y).all()
    d = out / ep
    d.mkdir()
    pred = np.zeros(750)
    info = []
    for name, weight in sorted(weights.items()):
        c = candidates[name]
        x = views[c["view"]]
        model = lgb.LGBMRegressor(**c["params"])
        model.fit(x[itr], y)
        p = model.predict(x[itest])
        assert np.isfinite(p).all()
        pred += weight * p
        model.booster_.save_model(str(d / f"{name}.txt"))
        info.append(dict(name=name, weight=weight, **c, n_train=len(y), n_features=x.shape[1]))
        print(
            json.dumps(dict(endpoint=ep, candidate=name, weight=weight, n_train=len(y))), flush=True
        )
    rows.append(
        testids.with_columns(
            pl.lit(ep).alias("endpoint"),
            pl.Series("prediction", pred),
            pl.lit("full_train").alias("split_id"),
            pl.lit(-1).alias("fold"),
        )
    )
    memberships.append(
        train.select("compound_id", "endpoint", "value").with_columns(
            pl.lit("full_train").alias("role")
        )
    )
    (d / "selected.json").write_text(json.dumps(info, indent=2))
    manifest.append(dict(endpoint=ep, n_train=len(y), selected=info))
pl.concat(rows).write_parquet(out / "test-predictions.parquet")
pl.concat(memberships).write_parquet(out / "membership.parquet")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            config=cfg,
            method=(
                "Average the five saved inner-selected weights per endpoint; "
                "refit each positively weighted original candidate on all endpoint labels. "
                "Fixed deployment seed 42; no test or new score-based selection."
            ),
            endpoints=manifest,
            feature_columns=feature_columns,
            n_train_labels=labels.height,
            n_test=testids.height,
            test_train_overlap=0,
        ),
        indent=2,
    )
)
