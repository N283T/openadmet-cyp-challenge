"""Independent statistical, split, checkpoint and numerical observation-model audit."""

import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
import torch.nn.functional as functional
from sklearn.model_selection import GroupShuffleSplit


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

    def path(name):
        found = list(inp.glob("*-" + name))
        assert len(found) == 1, (name, found)
        return found[0]

    def read(name):
        return pl.read_parquet(path(name))

    spec = importlib.util.spec_from_file_location("official", path("custom_scoring_functions.py"))
    scorer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scorer)
    oof, baseline = read("observation-oof.parquet"), read("focus_2d6_combined_50pct.parquet")
    key = ["compound_id", "endpoint", "fold", "split_id"]
    baseline = baseline.sort(key)
    assert baseline.height == 6384
    expected_methods = [arm + suffix for arm in cfg["arms"] for suffix in ["_direct", "_lgbm"]]
    assert set(oof["method"]) == set(expected_methods)
    metrics = []
    original_metrics = pl.read_csv(path("metrics.csv"))
    for method in expected_methods:
        table = oof.filter(pl.col("method") == method).sort(key)
        assert table.height == 6384 and table.select(key).equals(baseline.select(key))
        for col in ["y_true", "conf_low", "conf_high", "is_neighbor"]:
            assert np.array_equal(table[col].to_numpy(), baseline[col].to_numpy())
        assert np.isfinite(table["prediction"].to_numpy()).all()
        for (endpoint,), part in table.group_by("endpoint"):
            score = float(
                scorer.rae_soft_threshold_absolute_error(
                    part["y_true"].to_numpy(),
                    part["prediction"].to_numpy(),
                    y_true_lower=part["conf_low"].to_numpy(),
                    y_true_upper=part["conf_high"].to_numpy(),
                )
            )
            previous = original_metrics.filter(
                (pl.col("method") == method) & (pl.col("endpoint") == endpoint)
            )["ST_RAE"].item()
            assert np.isclose(score, previous, atol=1e-12, rtol=0)
            metrics.append(dict(method=method, endpoint=endpoint, n=part.height, ST_RAE=score))
    pl.DataFrame(metrics).write_csv(out / "official-metrics.csv")
    emb, scores = read("external-embeddings.parquet"), read("external-scores.parquet")
    base = emb.join(scores, on="compound_id", validate="1:1").sort("compound_id")
    ids = base["compound_id"].to_list()
    assert len(ids) == len(set(ids)) == 6896
    features = [c for c in base.columns if c.startswith(("emb_", "score_"))]
    x = base.select(features).to_numpy().astype(np.float32)
    assert x.shape == (6896, 647)
    lookup = {cid: index for index, cid in enumerate(ids)}
    kinds, enzymes = ["direct", "tdi", "emax", "logfc"], ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
    targets = read("observation-targets.parquet")
    assert targets.select("compound_id", "kind", "enzyme").n_unique() == targets.height
    y = np.full((6896, 16), np.nan, dtype=np.float32)
    lo, hi = y.copy(), y.copy()
    for row in targets.iter_rows(named=True):
        i = lookup[row["compound_id"]]
        j = 4 * kinds.index(row["kind"]) + enzymes.index(row["enzyme"])
        y[i, j] = row["value"]
        lo[i, j] = row["conf_low"] if row["conf_low"] is not None else row["value"]
        hi[i, j] = row["conf_high"] if row["conf_high"] is not None else row["value"]
    a = base.select("compound_id").join(
        read("parent-assignments.parquet"), on="compound_id", how="left", validate="1:1"
    )
    folds = a["fold"].fill_null(-1).to_numpy()
    anchors = a["is_anchor"].fill_null(False).to_numpy()
    groups = a["group"].cast(pl.String).fill_null("unassigned").to_numpy()
    available = (folds >= 0) & np.isfinite(y).any(axis=1)
    memberships, history = read("memberships.parquet"), read("training-history.parquet")
    assert memberships.height == 6896 * 5
    assert memberships.select("compound_id", "outer_fold").n_unique() == memberships.height
    summary = json.loads(path("summary.json").read_text())
    torch.set_num_threads(2)
    verification = []
    for fold in range(5):
        m = memberships.filter(pl.col("outer_fold") == fold).sort("compound_id")
        assert m["compound_id"].to_list() == ids
        train = available & ((folds != fold) | anchors)
        valid = available & (folds == fold) & ~anchors
        assert np.array_equal(m["outer_train"], train)
        assert np.array_equal(m["outer_validation"], valid)
        it, iv = m["inner_train"].to_numpy(), m["inner_validation"].to_numpy()
        assert np.array_equal(it | iv, train) and not (it & iv).any()
        assert not ((it | iv) & valid).any()
        assert not set(groups[it]) & set(groups[iv])
        assert it[available & anchors].all()
        tridx = np.flatnonzero(train)
        anchor_groups = set(groups[tridx[anchors[tridx]]])
        fixed = np.array([i for i in tridx if groups[i] in anchor_groups], dtype=int)
        nonanchor = np.array([i for i in tridx if groups[i] not in anchor_groups], dtype=int)
        ti, vi = next(
            GroupShuffleSplit(n_splits=1, test_size=0.1, random_state=42 + fold).split(
                nonanchor, groups=groups[nonanchor]
            )
        )
        assert set(np.flatnonzero(it)) == set(nonanchor[ti]) | set(fixed)
        assert set(np.flatnonzero(iv)) == set(nonanchor[vi])
        xm, xs = x[tridx].mean(0), np.maximum(x[tridx].std(0), 1e-6)
        mu, sd = np.nanmean(y[tridx], axis=0), np.maximum(np.nanstd(y[tridx], axis=0), 0.05)
        denominator = np.array(
            [
                (np.maximum(lo[tridx, j] - mu[j], 0) + np.maximum(mu[j] - hi[tridx, j], 0))[
                    np.isfinite(y[tridx, j])
                ].mean()
                for j in range(4)
            ],
            dtype=np.float32,
        )
        denominator = np.maximum(denominator, 0.05)
        for arm in cfg["arms"]:
            checkpoint = torch.load(
                path(f"{arm}-fold{fold}.pt"), map_location="cpu", weights_only=False
            )
            assert checkpoint["feature_columns"] == features and checkpoint["arm"] == arm
            for field, expected in [
                ("x_mean", xm),
                ("x_std", xs),
                ("head_mean", mu),
                ("head_std", sd),
                ("direct_denominator", denominator),
            ]:
                assert np.array_equal(checkpoint[field], expected), (
                    arm,
                    fold,
                    field,
                    np.max(np.abs(checkpoint[field] - expected)),
                )
            state = checkpoint["state_dict"]
            assert np.array_equal(state["mu"].numpy(), mu) and np.array_equal(
                state["sd"].numpy(), sd
            )
            h = history.filter((pl.col("arm") == arm) & (pl.col("fold") == fold))
            select = h.filter(pl.col("phase") == "select").sort("epoch")
            best, epoch = float("inf"), None
            for record in select.iter_rows(named=True):
                if record["inner_ST_RAE"] < best - 1e-6:
                    best, epoch = record["inner_ST_RAE"], record["epoch"]
            refit = h.filter(pl.col("phase") == "refit").sort("epoch")
            assert checkpoint["selected_epochs"] == epoch
            assert refit["epoch"].to_list() == list(range(1, epoch + 1))
            assert refit["inner_ST_RAE"].null_count() == refit.height
            fit = [r for r in summary["fits"] if r["arm"] == arm and r["fold"] == fold]
            assert len(fit) == 1 and fit[0]["epochs"] == epoch
            with torch.inference_mode():
                z = torch.tensor((x - xm) / xs)
                hidden = functional.relu(
                    functional.linear(z, state["net.0.weight"], state["net.0.bias"])
                )
                latent = functional.relu(
                    functional.linear(hidden, state["net.3.weight"], state["net.3.bias"])
                )
                prediction = (
                    functional.linear(latent, state["head.weight"], state["head.bias"])
                    * state["sd"]
                    + state["mu"]
                ).numpy()
            saved = read(f"{arm}-fold{fold}-embeddings.parquet").sort("compound_id")
            assert saved["compound_id"].to_list() == ids
            actual_latent = saved.select([f"emb_{j:03d}" for j in range(128)]).to_numpy()
            error = float(np.max(np.abs(latent.numpy() - actual_latent)))
            # Float32 CPU/GPU GEMMs can differ near zero after ReLU.
            assert np.allclose(latent.numpy(), actual_latent, atol=2e-4, rtol=2e-5), error
            direct = oof.filter((pl.col("method") == arm + "_direct") & (pl.col("fold") == fold))
            predicted = np.array(
                [
                    prediction[lookup[r["compound_id"]], enzymes.index(r["endpoint"].split("_")[0])]
                    for r in direct.iter_rows(named=True)
                ]
            )
            p_error = float(np.max(np.abs(predicted - direct["prediction"].to_numpy())))
            assert np.allclose(predicted, direct["prediction"], atol=2e-5, rtol=2e-5), p_error
            verification.append(
                dict(
                    arm=arm,
                    fold=fold,
                    all_train_statistics_exact=True,
                    selected_epochs=epoch,
                    refit_epochs=refit.height,
                    outer_train=int(train.sum()),
                    outer_validation=int(valid.sum()),
                    latent_max_abs_error=error,
                    prediction_max_abs_error=p_error,
                )
            )
    (out / "verification.json").write_text(
        json.dumps(
            dict(
                model_run=cfg["model_run"],
                methods=6,
                rows_per_method=6384,
                memberships=memberships.height,
                checks=verification,
            ),
            indent=2,
        )
    )
    print(
        json.dumps(dict(checkpoints_verified=len(verification), methods=6, oof_rows=oof.height)),
        flush=True,
    )


if __name__ == "__main__":
    main()
