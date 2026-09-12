"""Independent replay of greedy selections, thresholds and held-out predictions."""

import itertools
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import rankdata
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
audit = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
cfg = json.loads((inp / audit["config"]).read_text())


def path(name):
    exact = [p for p in inp.iterdir() if p.name.partition("-")[2] == name]
    if len(exact) == 1:
        return exact[0]
    candidates = [p for p in inp.iterdir() if p.name.endswith("-" + name)]
    assert len(candidates) == 1, (name, candidates)
    return candidates[0]


canonical = pl.read_parquet(inp / audit["canonical"])
keys = ["compound_id", "endpoint", "fold"]
banks = {}
for source in cfg["sources"]:
    for kind in ["oof", "inner", "test"]:
        frame = pl.read_parquet(path(source[kind]))
        if kind == "inner":
            if "fold" in frame.columns:
                frame = frame.drop("fold")
            frame = frame.rename({"outer_fold": "fold"})
        prediction = "probability" if "probability" in frame.columns else "prediction"
        for (method,), t in frame.group_by("method"):
            name = source["prefix"] + "::" + method
            columns = keys + [pl.col(prediction).cast(pl.Float64).alias("prediction")]
            if kind != "test":
                columns += [pl.col("y_true").cast(pl.Int64)]
            t = t.select(columns).sort(keys)
            assert t.height == (7500 if kind == "test" else 5081)
            assert not t.select(keys).is_duplicated().any()
            if kind != "test":
                joined = t.join(
                    canonical.select(
                        "compound_id",
                        "endpoint",
                        pl.col("fold").alias("actual_fold"),
                        pl.col("y_true").alias("actual_y"),
                    ),
                    on=["compound_id", "endpoint"],
                    validate="1:1",
                )
                assert joined.height == 5081
                assert (joined["y_true"] == joined["actual_y"]).all()
                expected = (joined["fold"] + 1) % 5 if kind == "inner" else joined["fold"]
                assert (joined["actual_fold"] == expected).all()
            banks[kind, name] = t
names = sorted(name for kind, name in banks if kind == "oof")
assert len(names) in [17, 18, 22, 24]
selections = json.loads(path("selections.json").read_text())
fitting = pl.read_csv(path("full-oof-fitting-only.csv"))
saved = {
    name: pl.read_parquet(path(file))
    for name, file in [
        ("oof", "oof.parquet"),
        ("test", "outer-weighted-fold-test.parquet"),
        ("deploy", "deployment.parquet"),
    ]
}
assert saved["oof"].height == len(cfg["recipes"]) * 5081
assert saved["test"].height == len(cfg["recipes"]) * 7500
assert saved["deploy"].height == len(cfg["recipes"]) * 1500


def cutscore(y, p):
    # Ascending cumulative counts, independently from producer descending sort.
    order = np.argsort(p, kind="stable")
    values = p[order]
    cuts, starts = np.unique(values, return_index=True)
    prior_positive = np.r_[0, np.cumsum(y[order])][starts]
    tp = y.sum() - prior_positive
    fp = len(y) - starts - tp
    fn = y.sum() - tp
    tn = len(y) - y.sum() - fp
    den = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    scores = np.divide(tp * tn - fp * fn, den, out=np.zeros_like(den, dtype=float), where=den > 0)
    options = list(zip(scores, -abs(cuts - 0.5), cuts, strict=True))
    options.append((matthews_corrcoef(y, p >= 0.5), 0.0, 0.5))
    best = max(options)
    return float(best[0]), float(best[2])


def objective(y, p, mode):
    if mode == "mcc":
        return cutscore(y, p)[0]
    if mode == "brier":
        return -np.mean((p - y) ** 2)
    ranks = rankdata(p)
    positive = y.sum()
    return (ranks[y == 1].sum() - positive * (positive + 1) / 2) / (positive * (len(y) - positive))


def bootstrap_sparse(x, y):
    counts = np.zeros(x.shape[1], dtype=int)
    running = np.zeros(len(y))
    best = -np.inf
    best_weights = None
    for step in range(1, cfg["steps"] + 1):
        allowed = [j for j in range(x.shape[1]) if counts[j] or np.count_nonzero(counts) < 3]
        scores = [objective(y, (running + x[:, j]) / step, "mcc") for j in allowed]
        chosen = allowed[np.argmax(scores)]
        counts[chosen] += 1
        running += x[:, chosen]
        if max(scores) > best + 1e-12:
            best = max(scores)
            best_weights = counts.copy() / step
    return best_weights


checks = []
for selected in selections:
    ep, recipe, fold = selected["endpoint"], selected["recipe"], selected["fold"]
    views = {}
    for kind in ["oof", "inner", "test"]:
        parts = [banks[kind, name].filter(pl.col("endpoint") == ep) for name in names]
        assert all(p.select(keys).equals(parts[0].select(keys)) for p in parts)
        views[kind] = (parts[0], np.column_stack([p["prediction"].to_numpy() for p in parts]))
    kind = "oof" if fold == -1 else "inner"
    ref, all_x = views[kind]
    mask = np.ones(ref.height, dtype=bool) if fold == -1 else ref["fold"].to_numpy() == fold
    x, y = all_x[mask], ref.filter(pl.Series(mask))["y_true"].to_numpy()
    weight = np.array([selected["weights"][name] for name in names])
    assert np.all(weight >= 0) and abs(weight.sum() - 1) < 1e-12
    best = -np.inf
    expected = None
    if recipe == "fixed_historical":
        expected = np.zeros(len(names))
        expected[names.index("history::historical_stack_mordred")] = 1
        best = objective(y, x @ expected, "mcc")
    elif recipe == "bagged_sparse3":
        rng = np.random.default_rng(cfg["bagging_seed"])
        batches = []
        assert len(selected["path"]) == cfg["bagging_samples"]
        for item in selected["path"]:
            sample = np.concatenate(
                [
                    rng.choice(
                        np.flatnonzero(y == value), size=int((y == value).sum()), replace=True
                    )
                    for value in [0, 1]
                ]
            )
            reconstructed = bootstrap_sparse(x[sample], y[sample])
            assert np.max(abs(reconstructed - np.asarray(item["weights"]))) < 1e-12
            sample_mcc = cutscore(y[sample], x[sample] @ reconstructed)[0]
            assert abs(sample_mcc - item["fit_MCC"]) < 1e-12
            batches.append(reconstructed)
        mean = np.mean(batches, axis=0)
        ranked = sorted(range(len(names)), key=lambda j: (-mean[j], names[j]))[:3]
        expected = np.zeros(len(names))
        expected[ranked] = mean[ranked]
        expected /= expected.sum()
        best = objective(y, x @ expected, "mcc")
    elif recipe.startswith("subset"):
        count = int(recipe[-1])
        for members in itertools.combinations(range(len(names)), count):
            value = objective(y, x[:, members].mean(axis=1), "mcc")
            if value > best + 1e-12:
                best = value
                expected = np.zeros(len(names))
                expected[list(members)] = 1 / count
    else:
        mode = recipe.split("_")[1]
        counts = np.zeros(len(names), dtype=int)
        running = np.zeros(len(y))
        cap = 3 if recipe.endswith("sparse3") else len(names)
        for step, item in enumerate(selected["path"], 1):
            allowed = [
                j for j in range(len(names)) if counts[j] > 0 or np.count_nonzero(counts) < cap
            ]
            scores = [objective(y, (running + x[:, j]) / step, mode) for j in allowed]
            chosen = allowed[np.argmax(scores)]
            assert item["added"] == names[chosen]
            assert abs(item["objective"] - max(scores)) < 1e-12
            counts[chosen] += 1
            running += x[:, chosen]
            if max(scores) > best + 1e-12:
                best = max(scores)
                expected = counts.copy() / step
        assert len(selected["path"]) == cfg["steps"]
    assert np.max(abs(expected - weight)) < 1e-12
    assert abs(best - selected["objective"]) < 1e-12
    score, cutoff = cutscore(y, x @ weight)
    assert abs(cutoff - selected["threshold"]) < 1e-12
    assert abs(score - selected["fit_MCC"]) < 1e-12
    if fold == -1:
        recorded = fitting.filter((pl.col("endpoint") == ep) & (pl.col("method") == recipe))
        assert recorded.height == 1
        p = x @ weight
        values = dict(
            MCC=matthews_corrcoef(y, p >= cutoff),
            AUROC=roc_auc_score(y, p),
            AP=average_precision_score(y, p),
            Brier=np.mean((y - p) ** 2),
            n=len(y),
            positives=int((p >= cutoff).sum()),
        )
        for key, value in values.items():
            assert abs(value - recorded[key][0]) < 1e-12
    kinds = ["deploy"] if fold == -1 else ["oof", "test"]
    maximum = 0.0
    for destination in kinds:
        raw_kind = "test" if destination == "deploy" else destination
        ref, matrix = views[raw_kind]
        expected_table = ref.select(keys).with_columns(pl.Series("expected", matrix @ weight))
        if destination == "deploy":
            expected_table = (
                expected_table.group_by("compound_id", "endpoint")
                .agg(pl.col("expected").mean())
                .with_columns(pl.lit(-1).alias("fold"))
            )
        else:
            expected_table = expected_table.filter(pl.col("fold") == fold)
        actual = (
            saved[destination]
            .filter(
                (pl.col("endpoint") == ep) & (pl.col("method") == recipe) & (pl.col("fold") == fold)
            )
            .sort(keys)
        )
        expected_table = expected_table.sort(keys)
        assert actual.select(keys).equals(expected_table.select(keys))
        error = np.max(abs(actual["prediction"].to_numpy() - expected_table["expected"].to_numpy()))
        assert error < 1e-12
        assert (abs(actual["threshold"] - cutoff) < 1e-12).all()
        maximum = max(maximum, float(error))
        if destination == "deploy":
            assert (actual["binary"] == (actual["prediction"] >= actual["threshold"])).all()
    checks.append(
        dict(
            endpoint=ep,
            recipe=recipe,
            fold=fold,
            max_error=maximum,
            members=int(np.count_nonzero(weight)),
        )
    )
    print(ep, recipe, fold, "passed", flush=True)
metrics = pl.read_csv(path("metrics.csv"))
for row in metrics.iter_rows(named=True):
    table = saved["oof"].filter(
        (pl.col("endpoint") == row["endpoint"]) & (pl.col("method") == row["method"])
    )
    if row["fold"] != -1:
        table = table.filter(pl.col("fold") == row["fold"])
    y, p = table["y_true"].to_numpy(), table["prediction"].to_numpy()
    binary = p >= table["threshold"].to_numpy()
    values = dict(
        MCC=matthews_corrcoef(y, binary),
        AUROC=roc_auc_score(y, p),
        AP=average_precision_score(y, p),
        Brier=np.mean((y - p) ** 2),
        n=len(y),
        positives=int(binary.sum()),
    )
    for key, value in values.items():
        assert abs(value - row[key]) < 1e-12
pl.DataFrame(checks).write_csv(out / "verification.csv")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            source_run=audit["source_run"],
            candidates=len(names),
            selections=len(checks),
            metrics=metrics.height,
            passed=True,
        ),
        indent=2,
    )
)
