"""Independently inspect TabICL serialized training arrays and canonical prediction alignment."""

import json
import os
import pickle
from pathlib import Path

import joblib
import numpy as np
import polars as pl
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
audit = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
cfg = json.loads((inp / audit["config"]).read_text())


def path(name):
    exact = [p for p in inp.iterdir() if p.name.partition("-")[2] == name]
    assert len(exact) == 1, (name, exact)
    return exact[0]


class StateOnly:
    def __setstate__(self, state):
        self.__dict__.update(state)


class StateLoader(pickle.Unpickler):
    def find_class(self, module, name):
        if name == "TabICLClassifier" and module.startswith("tabicl."):
            return StateOnly
        return super().find_class(module, name)


labels = pl.read_parquet(path(cfg["labels"]))
assignments = pl.read_parquet(path(cfg["assignments"]))
canonical = pl.read_parquet(inp / audit["canonical"]).filter(pl.col("strategy") == "umap")
assignment_columns = ["compound_id", "fold", "group"]
assert assignments.select(assignment_columns).sort("compound_id").equals(
    canonical.select(assignment_columns).sort("compound_id")
)
labels = labels.select("compound_id", "endpoint", "value").join(
    assignments.select("compound_id", "fold"), on="compound_id", validate="m:1"
)
oof = pl.read_parquet(path("outer-oof.parquet"))
inner = pl.read_parquet(path("inner-probabilities.parquet"))
test = pl.read_parquet(path("fold-test-probabilities.parquet"))
mean = pl.read_parquet(path("test-probabilities.parquet"))
assert oof.height == inner.height == len(cfg["views"]) * 5081
assert test.height == len(cfg["views"]) * 7500
assert mean.height == len(cfg["views"]) * 1500
checks = []
for view in cfg["views"]:
    feature = pl.concat([pl.read_parquet(path(name)) for name in view["files"]])
    columns = [name for name in feature.columns if name.startswith(view["prefix"])]
    for endpoint in cfg["expected_counts"]:
        frame = (
            labels.filter(pl.col("endpoint") == endpoint)
            .join(feature, on="compound_id", validate="1:1")
            .sort("compound_id")
        )
        raw = frame.select(columns).to_numpy()
        for fold in range(5):
            for phase in ["outer", "inner"]:
                name = f"{view['name']}-{endpoint}-fold{fold}-{phase}"
                member = pl.read_parquet(path(name + "-membership.parquet")).sort("compound_id")
                assert member.select("compound_id", "fold").equals(
                    frame.select("compound_id", "fold")
                )
                train = frame["fold"].to_numpy() != fold
                validation_fold = fold if phase == "outer" else (fold + 1) % 5
                if phase == "inner":
                    train &= frame["fold"].to_numpy() != validation_fold
                val = frame["fold"].to_numpy() == validation_fold
                np.testing.assert_array_equal(member["training"], train)
                np.testing.assert_array_equal(member["validation"], val)
                prepared = joblib.load(path(name + "-pca.joblib"))
                scaler, pca = prepared["scaler"], prepared["pca"]
                assert prepared["columns"] == columns
                np.testing.assert_allclose(
                    scaler.mean_, np.mean(raw[train], axis=0), atol=1e-6, rtol=2e-6
                )
                np.testing.assert_allclose(
                    scaler.var_, np.var(raw[train], axis=0), atol=1e-6, rtol=2e-6
                )
                assert scaler.n_samples_seen_ == int(train.sum())
                standardized = scaler.transform(raw[train])
                np.testing.assert_allclose(
                    pca.mean_, standardized.mean(axis=0), atol=1e-7, rtol=1e-5
                )
                assert pca.n_samples_ == int(train.sum()) and pca.n_components_ == 256
                projected = pca.transform(standardized).astype(np.float32)
                with path(name + ".tabicl_fit").open("rb") as stream:
                    state = StateLoader(stream).load()
                np.testing.assert_array_equal(state.classes_, [0, 1])
                assert state.n_estimators == 8 and state.random_state == 42 + fold
                assert state.allow_auto_download is False
                assert Path(state.model_path).name == cfg["weights"]
                stored_x = state.ensemble_generator_.X_
                np.testing.assert_allclose(
                    stored_x,
                    state.X_encoder_.transform(projected),
                    atol=1e-10,
                    rtol=1e-7,
                )
                np.testing.assert_array_equal(
                    state.ensemble_generator_.y_, frame["value"].to_numpy()[train]
                )
                selected = oof if phase == "outer" else inner
                key = "fold" if phase == "outer" else "outer_fold"
                selected = selected.filter(
                    (pl.col("method") == view["name"])
                    & (pl.col("endpoint") == endpoint)
                    & (pl.col(key) == fold)
                ).sort("compound_id")
                assert selected["compound_id"].equals(frame.filter(pl.Series(val))["compound_id"])
                np.testing.assert_array_equal(selected["y_true"], frame["value"].to_numpy()[val])
                assert selected["probability"].is_between(0, 1).all()
                checks.append(
                    dict(
                        method=view["name"],
                        endpoint=endpoint,
                        fold=fold,
                        phase=phase,
                        train=int(train.sum()),
                        validation=int(val.sum()),
                    )
                )
                print(name, "passed", flush=True)


def cutoff(y, p):
    options = [
        (matthews_corrcoef(y, p >= t), -abs(t - 0.5), -t, t) for t in np.linspace(0.01, 0.99, 99)
    ]
    return max(options)[3]


rows = []
for family, outer, inn in [
    ("tabicl", oof, inner),
    ("tabpfn", pl.read_parquet(inp / audit["pfn"]), pl.read_parquet(inp / audit["pfn_inner"])),
]:
    for (method, endpoint), t in outer.group_by("method", "endpoint"):
        pieces = []
        for fold in range(5):
            v = t.filter(pl.col("fold") == fold)
            s = inn.filter(
                (pl.col("method") == method)
                & (pl.col("endpoint") == endpoint)
                & (pl.col("outer_fold") == fold)
            )
            assert s.height and set(v["compound_id"]).isdisjoint(s["compound_id"])
            cut = cutoff(s["y_true"].to_numpy(), s["probability"].to_numpy())
            if family == "tabicl":
                np.testing.assert_allclose(v["inner_threshold"], cut, atol=1e-12, rtol=0)
            pieces.append(v.with_columns((pl.col("probability") >= cut).alias("audited_binary")))
        t = pl.concat(pieces)
        y, p = t["y_true"].to_numpy(), t["probability"].to_numpy()
        rows.append(
            dict(
                family=family,
                method=method,
                endpoint=endpoint,
                MCC05=matthews_corrcoef(y, p >= 0.5),
                MCC99inner=matthews_corrcoef(y, t["audited_binary"]),
                AUROC=roc_auc_score(y, p),
                AP=average_precision_score(y, p),
                n=len(y),
            )
        )
rebuilt = test.group_by("compound_id", "endpoint", "method").agg(
    pl.col("probability").mean().alias("rebuilt"), pl.len().alias("n")
)
paired = mean.join(rebuilt, on=["compound_id", "endpoint", "method"], validate="1:1")
assert (paired["n"] == 5).all()
np.testing.assert_allclose(paired["probability"], paired["rebuilt"], atol=1e-12, rtol=0)
pl.DataFrame(checks).write_csv(out / "fit-verification.csv")
pl.DataFrame(rows).sort("endpoint", "method").write_csv(out / "matched-metrics.csv")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            run=audit["run"],
            fits=len(checks),
            canonical_assignment_rows=canonical.height,
            canonical_assignments_exact=True,
            passed=True,
        ),
        indent=2,
    )
)
