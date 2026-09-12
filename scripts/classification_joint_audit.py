"""Audit multitask TDI exclusions, target scalers, epoch lineage and keyed outputs."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
from sklearn.metrics import matthews_corrcoef

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(key):
    return pl.read_parquet(inp / cfg[key])


labels, aux, assigned, compounds = (
    read("labels"),
    read("aux"),
    read("assignments"),
    read("compounds"),
)
assigned = assigned.filter(pl.col("strategy") == "umap")
class_heads = ["CYP2D6_is_TDI", "CYP3A4_is_TDI"]
aux = aux.filter(
    (pl.col("kind") == "logfc")
    | ((pl.col("kind") == "tdi") & pl.col("enzyme").is_in(["CYP2D6", "CYP3A4"]))
).select("compound_id", (pl.col("enzyme") + "_" + pl.col("kind")).alias("endpoint"), "value")
direct = labels.filter(pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition")).select(
    "compound_id", "endpoint", "value"
)
aux = pl.concat([aux, direct])
heads = class_heads + sorted(aux["endpoint"].unique())
truth = pl.concat(
    [
        labels.filter(pl.col("endpoint").is_in(class_heads)).select(
            "compound_id", "endpoint", "value"
        ),
        aux,
    ]
)
table = (
    compounds.select("compound_id")
    .sort("compound_id")
    .join(
        assigned.select("compound_id", "fold"),
        on="compound_id",
        how="left",
        validate="1:1",
        maintain_order="left",
    )
    .join(
        truth.pivot(on="endpoint", index="compound_id", values="value"),
        on="compound_id",
        how="left",
        validate="1:1",
        maintain_order="left",
    )
)
assert table.height == 6896 and len(heads) == 12
values = table.select(heads).to_numpy().astype(np.float32)
folds = table["fold"].fill_null(-1).to_numpy()
contract = json.loads((inp / cfg["contract"]).read_text())
assert contract["heads"] == heads
assert contract["counts"] == np.isfinite(values).sum(0).tolist()
oof, inner, test, foldtest = read("oof"), read("inner"), read("test"), read("foldtest")
assert oof.height == inner.height == 10162 and test.height == 3000 and foldtest.height == 15000
assert not foldtest["y_true"].is_finite().any()
assert oof.unique(["compound_id", "endpoint", "method"]).height == 10162
assert test.unique(["compound_id", "endpoint", "method"]).height == 3000
audit_rows = []
for item in cfg["fits"]:
    recipe, outer = item["recipe"], item["fold"]
    nheads = 2 if recipe == "tdi_only" else 12
    innerfold = (outer + 1) % 5
    checkpoints = {}
    for phase in ["inner", "refit"]:
        expected = (folds >= 0) & (folds != outer) & np.isfinite(values[:, :nheads]).any(1)
        if phase == "inner":
            expected &= folds != innerfold
        membership = pl.read_parquet(inp / item[f"{phase}_membership"]).sort("compound_id")
        np.testing.assert_array_equal(membership["compound_id"], table["compound_id"])
        np.testing.assert_array_equal(membership["training"], expected)
        np.testing.assert_array_equal(membership["outer_holdout"], folds == outer)
        np.testing.assert_array_equal(
            membership["inner_validation"], (folds == innerfold) & (phase == "inner")
        )
        scaling = json.loads((inp / item[f"{phase}_scaling"]).read_text())
        mean, std = np.zeros(nheads, dtype=np.float32), np.ones(nheads, dtype=np.float32)
        if nheads > 2:
            mean[2:] = np.nanmean(values[expected, 2:nheads], axis=0)
            std[2:] = np.maximum(np.nanstd(values[expected, 2:nheads], axis=0), 1e-6)
        np.testing.assert_array_equal(scaling["mean"], mean)
        np.testing.assert_array_equal(scaling["std"], std)
        assert scaling["train_counts"] == np.isfinite(values[expected, :nheads]).sum(0).tolist()
        checkpoint = torch.load(
            inp / item[f"{phase}_checkpoint"], map_location="cpu", weights_only=False
        )
        checkpoints[phase] = {k: checkpoint[k] for k in checkpoint if k != "state_dict"}
        del checkpoint
    history = json.loads((inp / item["inner_history"]).read_text())
    best, selected = float("inf"), None
    for record in history:
        if record["inner_BCE"] < best - 1e-6:
            best, selected = record["inner_BCE"], record["epoch"]
    refit_history = json.loads((inp / item["refit_history"]).read_text())
    assert (
        selected
        == checkpoints["inner"]["epoch"]
        == checkpoints["refit"]["epochs"]
        == len(refit_history)
    )
    assert [row["epoch"] for row in refit_history] == list(range(1, selected + 1))
    for j, endpoint in enumerate(class_heads):
        selected_inner = inner.filter(
            (pl.col("method") == recipe)
            & (pl.col("endpoint") == endpoint)
            & (pl.col("outer_fold") == outer)
        )
        canonical = labels.filter(pl.col("endpoint") == endpoint).join(
            assigned.select("compound_id", "fold"), on="compound_id", validate="1:1"
        )
        assert set(selected_inner["compound_id"]) == set(
            canonical.filter(pl.col("fold") == innerfold)["compound_id"]
        )
        candidates = np.unique(np.r_[0.5, np.linspace(0.05, 0.95, 91)])
        scores = [
            matthews_corrcoef(selected_inner["y_true"], selected_inner["probability"] >= t)
            for t in candidates
        ]
        threshold = float(candidates[np.argmax(scores)])
        assert threshold == checkpoints["inner"]["thresholds"][j]
        predictions = oof.filter(
            (pl.col("method") == recipe)
            & (pl.col("endpoint") == endpoint)
            & (pl.col("fold") == outer)
        )
        assert set(predictions["compound_id"]) == set(
            canonical.filter(pl.col("fold") == outer)["compound_id"]
        )
        compare = predictions.join(
            canonical.select("compound_id", "value"), on="compound_id", validate="1:1"
        )
        np.testing.assert_array_equal(compare["y_true"], compare["value"])
        np.testing.assert_array_equal(
            predictions["threshold"], np.full(predictions.height, threshold)
        )
        np.testing.assert_array_equal(
            predictions["prediction"], (predictions["probability"] >= threshold).cast(pl.Int8)
        )
        audit_rows.append(
            dict(
                recipe=recipe,
                outer_fold=outer,
                endpoint=endpoint,
                selected_epoch=selected,
                threshold=threshold,
                oof_rows=predictions.height,
                inner_rows=selected_inner.height,
            )
        )
for recipe in ["tdi_only", "tdi_aux"]:
    assert oof.filter(pl.col("method") == recipe).height == 5081
    assert test.filter(pl.col("method") == recipe).height == 1500
    for endpoint in class_heads:
        subset = test.filter((pl.col("method") == recipe) & (pl.col("endpoint") == endpoint))
        assert set(subset["compound_id"]) == set(read("testmap")["compound_id"])
rebuilt = foldtest.group_by("compound_id", "endpoint", "method").agg(
    pl.col("probability").mean().alias("p"),
    pl.col("threshold").mean().alias("t"),
    pl.len().alias("count"),
)
paired = test.join(rebuilt, on=["compound_id", "endpoint", "method"], validate="1:1")
assert (paired["count"] == 5).all()
np.testing.assert_allclose(paired["probability"], paired["p"], rtol=0, atol=1e-7)
np.testing.assert_allclose(paired["threshold"], paired["t"], rtol=0, atol=1e-12)
assert oof["probability"].is_between(0, 1).all() and test["probability"].is_between(0, 1).all()
pl.DataFrame(audit_rows).write_csv(out / "fit-audit.csv")
(out / "summary.json").write_text(
    json.dumps(
        dict(
            model_run=cfg["model_run"],
            fits=20,
            outer_auxiliary_exclusion_exact=True,
            inner_exclusion_exact=True,
            scalers_exact=True,
            epoch_lineage_exact=True,
            thresholds_reconstructed=True,
            oof_rows=10162,
            test_rows=3000,
            test_mean_max_error=float((paired["probability"] - paired["p"]).abs().max()),
        ),
        indent=2,
    )
)
