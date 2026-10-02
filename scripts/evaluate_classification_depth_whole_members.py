"""Reconstruct submitted TDI OOF and screen whole CheMeleon depth members."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
from scipy.stats import rankdata
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def fold_rank(values, folds):
    result = np.empty(len(values), dtype=float)
    for fold in np.unique(folds):
        mask = folds == fold
        result[mask] = rankdata(values[mask], method="average") / (mask.sum() + 1)
    return result


def aligned(ids, frame, method, endpoint, column="prediction"):
    part = frame.filter((pl.col("endpoint") == endpoint) & (pl.col("method") == method))
    if "strategy" in part.columns:
        part = part.filter(pl.col("strategy") == "umap")
    result = ids.join(
        part.select("compound_id", pl.col(column).alias("value")),
        on="compound_id",
        validate="1:1",
        maintain_order="left",
    )
    assert result.height == ids.height and result["value"].is_finite().all(), method
    return result["value"].to_numpy()


def metrics(y, score, calls, prior):
    positive = y == 1
    tpr = float(calls[positive].mean())
    fpr = float(calls[~positive].mean())
    tp, fn = prior * tpr, prior * (1 - tpr)
    fp, tn = (1 - prior) * fpr, (1 - prior) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return {
        "AUROC": float(roc_auc_score(y, score)),
        "AP": float(average_precision_score(y, score)),
        "MCC_observed": float(matthews_corrcoef(y, calls)),
        "MCC_target_prior": float((tp * tn - fp * fn) / denominator) if denominator else 0.0,
    }


def fixed_calls(score, folds, rate):
    calls = np.zeros(len(score), dtype=bool)
    for fold in np.unique(folds):
        rows = np.flatnonzero(folds == fold)
        count = max(1, int(round(rate * len(rows))))
        calls[rows[np.argsort(-score[rows], kind="stable")[:count]]] = True
    return calls


depth_frame = pl.read_parquet(inp / "depth-oof.parquet").filter(pl.col("track") == "classification")
rows, prediction_frames = [], []

# D6: preserve the submitted fold-selected core, its fixed-model tail and the
# PCA/raw block. A depth readout enters the core as one complete member.
d6_saved = pl.read_parquet(inp / "d6-current.parquet").sort("compound_id")
d6_ids = d6_saved.select("compound_id")
d6_folds = d6_saved["fold"].to_numpy()
d6_y = d6_saved["y_true"].to_numpy().astype(int)
d6_endpoint = "CYP2D6_is_TDI"
d6_specs = {
    "lgbm::morgan_plain": ("d6-lgbm.parquet", "morgan_plain", "prediction"),
    "lgbm::mordred_plain": ("d6-lgbm.parquet", "mordred_plain", "prediction"),
    "lgbm::fp_mordred_plain": ("d6-lgbm.parquet", "fp_mordred_plain", "prediction"),
    "lgbm::external_class_plain": ("d6-lgbm.parquet", "external_class_plain", "prediction"),
    "lgbm::external_mixed_plain": ("d6-lgbm.parquet", "external_mixed_plain", "prediction"),
    "lgbm::external_class_mordred_plain": (
        "d6-lgbm.parquet",
        "external_class_mordred_plain",
        "prediction",
    ),
    "pfn::external_class_pca256": ("d6-pfn.parquet", "external_class_pca256", "probability"),
    "pfn::external_mixed_pca256": ("d6-pfn.parquet", "external_mixed_pca256", "probability"),
    "pfn::chemeleon_pca256": ("d6-pfn.parquet", "chemeleon_pca256", "probability"),
    "graph::tdi_only": ("d6-graph.parquet", "tdi_only", "probability"),
    "graph::tdi_aux": ("d6-graph.parquet", "tdi_aux", "probability"),
    "pfn::external_mixed_raw600": ("d6-raw600.parquet", "external_mixed_raw600", "probability"),
}
d6_atoms = {
    name: aligned(d6_ids, pl.read_parquet(inp / file), method, d6_endpoint, column)
    for name, (file, method, column) in d6_specs.items()
}
d6_fixed = pl.read_parquet(inp / "d6-fixed.parquet")
d6_tail = {
    label: fold_rank(aligned(d6_ids, d6_fixed, method, d6_endpoint), d6_folds)
    for label, method in {
        "logistic": "logistic_c0.05",
        "rbf": "rbf_svm_c1",
        "complement": "complement_nb_a1",
    }.items()
}
d6_choices = json.loads(
    pl.read_parquet(inp / "d6-selections.parquet").filter(
        pl.col("model") == "exhaustive_mcc::components_10"
    )["selection"][0]
)
d6_choices = [item for item in d6_choices if item["held_fold"] != "full"]


def d6_predict(depth=None, replace=False):
    result = np.empty(len(d6_y), dtype=float)
    for choice in d6_choices:
        members = choice["members"]
        core_members = []
        for member in members:
            if member == "pfn::external_mixed_pca256":
                block = (d6_atoms[member] + d6_atoms["pfn::external_mixed_raw600"]) / 2
                core_members.append(depth if replace else block)
            else:
                core_members.append(d6_atoms[member])
        if depth is not None and not replace:
            core_members.append(depth)
        core = fold_rank(np.mean(core_members, axis=0), d6_folds)
        live = (2 * core + d6_tail["logistic"]) / 3
        recipes = {
            "core": core,
            "core_logistic": live,
            "tail_rbf": (2 * live + d6_tail["rbf"]) / 3,
            "tail_complement": (2 * live + d6_tail["complement"]) / 3,
        }
        recipes["tail_average"] = (recipes["tail_rbf"] + recipes["tail_complement"]) / 2
        result[d6_folds == int(choice["held_fold"])] = recipes[choice["recipe"]][
            d6_folds == int(choice["held_fold"])
        ]
    return result


def depth_values(endpoint, ids, folds, levels):
    values = []
    for level in levels:
        part = depth_frame.filter(
            (pl.col("endpoint") == endpoint) & (pl.col("depth") == level)
        ).select("compound_id", "fold", "prediction")
        aligned_part = ids.join(part, on="compound_id", validate="1:1", maintain_order="left")
        assert np.array_equal(aligned_part["fold"].to_numpy(), folds)
        values.append(aligned_part["prediction"].to_numpy())
    return fold_rank(np.mean(values, axis=0), folds)


def backend_values(endpoint, ids, folds, filename):
    frame = pl.read_parquet(inp / filename).filter(pl.col("endpoint") == endpoint)
    part = ids.join(
        frame.select("compound_id", "fold", "probability"),
        on="compound_id",
        validate="1:1",
        maintain_order="left",
    )
    assert part.height == ids.height and np.array_equal(part["fold"].to_numpy(), folds)
    assert part["probability"].is_finite().all()
    return fold_rank(part["probability"].to_numpy(), folds)


d6_current = d6_predict()
assert np.max(np.abs(d6_current - d6_saved["prediction"].to_numpy())) < 1e-12
d6_scenarios = {"current": d6_current}
for levels in cfg["depth_sets"]:
    name = "+".join(map(str, levels))
    depth = depth_values(d6_endpoint, d6_ids, d6_folds, levels)
    d6_scenarios[f"add_depth_{name}"] = d6_predict(depth)
    d6_scenarios[f"replace_pca_raw_block_depth_{name}"] = d6_predict(depth, replace=True)
for backend, filename in cfg.get("backends", {}).items():
    depth = backend_values(d6_endpoint, d6_ids, d6_folds, filename)
    d6_scenarios[f"add_d1_{backend}"] = d6_predict(depth)
    d6_scenarios[f"replace_pca_raw_block_d1_{backend}"] = d6_predict(depth, replace=True)
if cfg.get("trunk_predictions"):
    trunk = pl.read_parquet(inp / cfg["trunk_predictions"]["d6"])
    for name in sorted(trunk["model"].unique()):
        part = trunk.filter((pl.col("model") == name) & (pl.col("role") == "train"))
        part = d6_saved.select("compound_id", "fold").join(
            part.select("compound_id", "fold", "prediction"),
            on=["compound_id", "fold"],
            validate="1:1",
            maintain_order="left",
        )
        assert part.height == d6_saved.height and part["prediction"].is_finite().all()
        candidate = fold_rank(part["prediction"].to_numpy(), d6_folds)
        d6_scenarios[f"standalone::{name}"] = candidate
        if name.startswith("molecular_only"):
            continue
        d6_scenarios[f"add::{name}"] = d6_predict(candidate)
        d6_scenarios[f"replace_pca_raw::{name}"] = d6_predict(candidate, replace=True)

# A4: reconstruct the seven rank-normalized atomic members and the saved
# held-fold-selected seven-member OOF before addition/replacement.
a4_saved = (
    pl.read_parquet(inp / "a4-oof.parquet")
    .filter(pl.col("model") == "exhaustive_mcc::components_7")
    .sort("compound_id")
)
a4_ids = a4_saved.select("compound_id")
a4_folds = a4_saved["fold"].to_numpy()
a4_y = a4_saved["y_true"].to_numpy().astype(int)
a4_endpoint = "CYP3A4_is_TDI"
a4_member = pl.read_parquet(inp / "a4-member.parquet").sort("compound_id")
assert a4_member["compound_id"].equals(a4_saved["compound_id"])
a4_atoms = {
    name: fold_rank(a4_member[column].to_numpy(), a4_folds)
    for name, column in {
        "historical": "history::historical_stack_mordred",
        "jeremy_tabicl_pca256": "jeremy8::tabicl_jeremy_pca256",
        "jeremy_tabpfn_pca256": "jeremy8::tabpfn_jeremy_pca256",
        "external_class_mordred": "lgbm::external_class_mordred_balanced",
    }.items()
}
raw = pl.read_parquet(inp / "a4-raw.parquet")
for name, method in {
    "jeremy_tabicl_raw600": "tabicl_jeremy_raw600",
    "jeremy_tabpfn_raw600": "tabpfn_jeremy_raw600",
}.items():
    a4_atoms[name] = fold_rank(aligned(a4_ids, raw, method, a4_endpoint, "probability"), a4_folds)
a4_atoms["pxr_admetica"] = fold_rank(
    aligned(
        a4_ids,
        pl.read_parquet(inp / "a4-aux.parquet"),
        "morgan_pxr_admetica_lgbm",
        a4_endpoint,
    ),
    a4_folds,
)
a4_selections = pl.read_parquet(inp / "a4-selections.parquet")
a4_choices = json.loads(
    a4_selections.filter(pl.col("model") == "exhaustive_mcc::components_7")["members"][0]
)
a4_choices = [item for item in a4_choices if item["held_fold"] != -1]


def a4_predict(depth=None, replaced=None):
    result = np.empty(len(a4_y), dtype=float)
    for choice in a4_choices:
        names = choice["members"]
        values = [depth if name == replaced else a4_atoms[name] for name in names]
        if depth is not None and replaced is None:
            values.append(depth)
        prediction = np.mean(values, axis=0)
        result[a4_folds == choice["held_fold"]] = prediction[a4_folds == choice["held_fold"]]
    return result


a4_current = a4_predict()
a4_error = np.abs(a4_current - a4_saved["prediction"].to_numpy())
assert np.max(a4_error) < 1e-8, (float(np.max(a4_error)), float(np.mean(a4_error)))
a4_scenarios = {"current": a4_current}
for levels in cfg["depth_sets"]:
    name = "+".join(map(str, levels))
    depth = depth_values(a4_endpoint, a4_ids, a4_folds, levels)
    a4_scenarios[f"add_depth_{name}"] = a4_predict(depth)
    a4_scenarios[f"replace_historical_depth_{name}"] = a4_predict(depth, "historical")
for backend, filename in cfg.get("backends", {}).items():
    depth = backend_values(a4_endpoint, a4_ids, a4_folds, filename)
    a4_scenarios[f"add_d1_{backend}"] = a4_predict(depth)
    a4_scenarios[f"replace_historical_d1_{backend}"] = a4_predict(depth, "historical")
if cfg.get("trunk_predictions"):
    trunk = pl.read_parquet(inp / cfg["trunk_predictions"]["a4"])
    for name in sorted(trunk["model"].unique()):
        part = trunk.filter((pl.col("model") == name) & (pl.col("role") == "train"))
        part = a4_saved.select("compound_id", "fold").join(
            part.select("compound_id", "fold", "prediction"),
            on=["compound_id", "fold"],
            validate="1:1",
            maintain_order="left",
        )
        assert part.height == a4_saved.height and part["prediction"].is_finite().all()
        candidate = fold_rank(part["prediction"].to_numpy(), a4_folds)
        a4_scenarios[f"standalone::{name}"] = candidate
        if name.startswith("molecular_only"):
            continue
        a4_scenarios[f"add::{name}"] = a4_predict(candidate)
        for replaced in a4_atoms:
            a4_scenarios[f"replace_{replaced}::{name}"] = a4_predict(candidate, replaced)

for target, ids, folds, y, scenarios in [
    ("d6", d6_ids, d6_folds, d6_y, d6_scenarios),
    ("a4", a4_ids, a4_folds, a4_y, a4_scenarios),
]:
    setting = cfg[target]
    if cfg.get("trunk_predictions"):
        selected = np.empty(len(y))
        selections = []
        candidates = {
            name: score for name, score in scenarios.items() if not name.startswith("standalone::")
        }
        for held in range(5):
            train = folds != held
            name = min(
                candidates,
                key=lambda name: (
                    -average_precision_score(y[train], candidates[name][train]),
                    name,
                ),
            )
            selected[~train] = candidates[name][~train]
            selections.append({"target": target, "held_fold": held, "scenario": name})
        scenarios["outer_selected_AP_placement"] = selected
        pl.DataFrame(selections).write_parquet(out / f"{target}-trunk-selections.parquet")
    for name, score in scenarios.items():
        calls = fixed_calls(score, folds, setting["call_rate"])
        rows.append(
            {
                "target": target,
                "scenario": name,
                "n": len(y),
                "calls": int(calls.sum()),
                **metrics(y, score, calls, setting["target_prior"]),
            }
        )
        prediction_frames.append(
            ids.with_columns(
                pl.lit(target).alias("target"),
                pl.Series("fold", folds),
                pl.Series("y_true", y),
                pl.lit(name).alias("scenario"),
                pl.Series("prediction", score),
                pl.Series("call", calls),
            )
        )
pl.DataFrame(rows).write_csv(out / "metrics.csv")
pl.concat(prediction_frames).write_parquet(out / "oof.parquet")
(out / "summary.json").write_text(json.dumps({"config": cfg}, indent=2) + "\n")
