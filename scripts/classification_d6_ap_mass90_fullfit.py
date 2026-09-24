"""Compose full-fit blind predictions for the unrestricted AP mass90 D6 pool."""

import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import rankdata
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
pool_cfg = json.loads((inp / cfg["pool_definition"]).read_text())


def read(key: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[key])


def unit_rank(values: np.ndarray) -> np.ndarray:
    return rankdata(values, method="average") / (len(values) + 1.0)


def fold_rank(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    result = np.empty(len(values), dtype=float)
    for fold in np.unique(folds):
        mask = folds == fold
        result[mask] = unit_rank(values[mask])
    return result


def fixed_calls(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    calls = np.zeros(len(values), dtype=bool)
    for fold in np.unique(folds):
        rows = np.flatnonzero(folds == fold)
        count = max(1, int(round(cfg["call_rate"] * len(rows))))
        calls[rows[np.argsort(-values[rows], kind="stable")[:count]]] = True
    return calls


def top_calls(values: np.ndarray, count: int) -> np.ndarray:
    calls = np.zeros(len(values), dtype=bool)
    calls[np.argsort(-values, kind="stable")[:count]] = True
    return calls


def fold_test_values(
    table: pl.DataFrame,
    test_ids: pl.DataFrame,
    predicate: pl.Expr,
    value_column: str = "prediction",
) -> np.ndarray:
    selected = table.filter(predicate)
    unique_folds = sorted(selected["fold"].unique().to_list())
    if unique_folds == [-1]:
        one = selected.select("compound_id", pl.col(value_column).alias("value"))
        joined = test_ids.join(one, on="compound_id", validate="1:1", maintain_order="left")
        if joined["value"].null_count():
            raise ValueError("Incomplete averaged test alignment")
        return unit_rank(joined["value"].to_numpy())
    fold_values = []
    for fold in unique_folds:
        one = selected.filter(pl.col("fold") == fold).select(
            "compound_id", pl.col(value_column).alias("value")
        )
        joined = test_ids.join(one, on="compound_id", validate="1:1", maintain_order="left")
        if joined["value"].null_count():
            raise ValueError("Incomplete test alignment")
        fold_values.append(unit_rank(joined["value"].to_numpy()))
    if len(fold_values) != 5:
        raise ValueError(f"Expected five test folds, found {len(fold_values)}")
    return unit_rank(np.mean(fold_values, axis=0))


legacy = read("legacy_bank").filter(pl.col("endpoint") == cfg["endpoint"])
reference = legacy.filter(pl.col("model") == "current_fixed_ratio").sort("compound_id")
ids = reference.select("compound_id")
folds = reference["fold"].to_numpy().astype(int)
y = reference["y_true"].to_numpy().astype(int)

outer = read("outer_selections").filter(pl.col("recipe") == cfg["selection_recipe"])
all_names = sorted(
    name for name in legacy["model"].unique().to_list() if name != "current_fixed_ratio"
) + [item["name"] for item in pool_cfg["jeremy_chemeleon_members"]]
weight_rows = []
for row in outer.iter_rows(named=True):
    weights = json.loads(row["weights"])
    weight_rows.append([weights.get(name, 0.0) for name in all_names])
mean_weight = np.mean(weight_rows, axis=0)
order = np.argsort(-mean_weight, kind="stable")
positive = order[mean_weight[order] > 0]
cumulative = np.cumsum(mean_weight[positive])
pool_count = int(np.searchsorted(cumulative, cfg["selection_weight_mass"], side="left") + 1)
pool_indices = positive[:pool_count]
pool_names = [all_names[index] for index in pool_indices]

oof_sources = {
    "lgbm": read("lgbm_oof"),
    "readout": read("readout_oof"),
    "tabpfn": read("tabpfn_oof"),
}
new_specs = {item["name"]: item for item in pool_cfg["jeremy_chemeleon_members"]}
oof_values: dict[str, np.ndarray] = {}
for name in pool_names:
    if name in new_specs:
        spec = new_specs[name]
        table = oof_sources[spec["source"]].filter(
            (pl.col("endpoint") == cfg["endpoint"])
            & (pl.col("split_id") == cfg["split_id"])
            & (pl.col("method") == spec["method"])
        )
    else:
        table = legacy.filter(pl.col("model") == name)
    joined = ids.join(
        table.select("compound_id", "prediction"),
        on="compound_id", validate="1:1", maintain_order="left",
    )
    if joined["prediction"].null_count():
        raise ValueError(f"Incomplete OOF alignment for {name}")
    oof_values[name] = fold_rank(joined["prediction"].to_numpy(), folds)

control_test = read("control_test")
control_test_one = control_test.filter(pl.col("model") == cfg["control_model"]).sort("compound_id")
test_ids = control_test_one.select("compound_id")
test_sources = {
    "lgbm": read("lgbm_test"),
    "readout": read("readout_test"),
    "tabpfn": read("tabpfn_test"),
}
test_values: dict[str, np.ndarray] = {}
for name in pool_names:
    if name in new_specs:
        spec = new_specs[name]
        test_values[name] = fold_test_values(
            test_sources[spec["source"]], test_ids,
            (pl.col("endpoint") == cfg["endpoint"])
            & (pl.col("split_id") == cfg["split_id"])
            & (pl.col("method") == spec["method"]),
        )
    elif name == "new::mace":
        test_values[name] = fold_test_values(
            read("mace_test"), test_ids,
            (pl.col("endpoint") == cfg["endpoint"])
            & (pl.col("strategy") == cfg["split_id"])
            & (pl.col("method") == "morgan_mace_off_emb_lr"),
        )
    elif name == "new::product_activity":
        test_values[name] = fold_test_values(
            read("activity_test"), test_ids,
            (pl.col("endpoint") == cfg["endpoint"])
            & (pl.col("strategy") == cfg["split_id"])
            & (pl.col("method") == "morgan_activity_endpoint_lgbm"),
        )
    elif name == "new::smartcyp":
        test_values[name] = fold_test_values(
            read("mechanism_test"), test_ids,
            (pl.col("endpoint") == cfg["endpoint"])
            & (pl.col("strategy") == cfg["split_id"])
            & (pl.col("method") == "morgan_smartcyp_lgbm"),
        )
    elif name == "pfn::external_mixed_pca256":
        test_values[name] = fold_test_values(
            read("pfn_test"), test_ids,
            (pl.col("endpoint") == cfg["endpoint"])
            & (pl.col("method") == "external_mixed_pca256"),
            "probability",
        )
    elif name == "complement_nb":
        compounds = read("compounds")
        mapping = read("test_mapping").select("compound_id").unique()
        ordered = test_ids.join(mapping, on="compound_id", validate="1:1").join(
            compounds.select("compound_id", "standardized_smiles"),
            on="compound_id", validate="1:1", maintain_order="left",
        )
        generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        features = np.stack([
            generator.GetFingerprintAsNumPy(Chem.MolFromSmiles(smiles))
            for smiles in ordered["standardized_smiles"]
        ]).astype(np.float32)
        fold_values = []
        for model_name in cfg["complement_models"]:
            model = joblib.load(inp / model_name)
            fold_values.append(unit_rank(model.predict_proba(features)[:, 1]))
        test_values[name] = unit_rank(np.mean(fold_values, axis=0))
    else:
        raise ValueError(f"No blind prediction source for {name}")

x_oof = np.column_stack([oof_values[name] for name in pool_names])
x_test = np.column_stack([test_values[name] for name in pool_names])
counts = np.zeros(len(pool_names), dtype=int)
total = np.zeros(len(y), dtype=float)
best_ap = -np.inf
best_weights = None
path = []
for step in range(1, cfg["steps"] + 1):
    values = []
    for index in range(len(pool_names)):
        prediction = (total + x_oof[:, index]) / step
        values.append(float(average_precision_score(y, prediction)))
    choice = int(np.argmax(values))
    total += x_oof[:, choice]
    counts[choice] += 1
    path.append({"step": step, "added": pool_names[choice], "AP": values[choice]})
    if values[choice] > best_ap + 1e-12:
        best_ap = values[choice]
        best_weights = counts.astype(float) / step
if best_weights is None:
    raise RuntimeError("Full-fit Caruana failed")

selected_indices = np.flatnonzero(best_weights)
selected_names = [pool_names[index] for index in selected_indices]
candidate_test_score = x_test @ best_weights
candidate_test_call = top_calls(candidate_test_score, cfg["test_positive_calls"])
adaptive = (
    read("adaptive_oof")
    .filter(pl.col("method") == cfg["adaptive_method"])
    .sort("compound_id")
)
control_oof = (
    read("control_oof")
    .filter(pl.col("model") == cfg["control_model"])
    .sort("compound_id")
)

decision_oof = pl.concat([
    control_oof.select("compound_id", "y_true", "model", "score", "call"),
    adaptive.select(
        "compound_id", "y_true",
        pl.lit(cfg["candidate_model"]).alias("model"),
        pl.col("prediction").alias("score"), "call",
    ),
])
decision_test = pl.concat([
    control_test_one.select("compound_id", "model", "score", "call"),
    test_ids.with_columns(
        pl.lit(cfg["candidate_model"]).alias("model"),
        pl.Series("score", candidate_test_score),
        pl.Series("call", candidate_test_call),
    ),
])
member_oof = pl.concat([
    ids.with_columns(
        pl.lit(name).alias("model"),
        pl.Series("prediction", oof_values[name]),
        pl.Series("call", fixed_calls(oof_values[name], folds)),
    )
    for name in selected_names
])
member_test = pl.concat([
    test_ids.with_columns(
        pl.lit(name).alias("model"),
        pl.Series("prediction", test_values[name]),
        pl.Series("call", top_calls(test_values[name], cfg["test_positive_calls"])),
    )
    for name in selected_names
])
decision_oof.write_parquet(out / "decision-oof.parquet")
decision_test.write_parquet(out / "decision-test.parquet")
member_oof.write_parquet(out / "member-oof.parquet")
member_test.write_parquet(out / "member-test.parquet")

honest_score = adaptive["prediction"].to_numpy()
honest_call = adaptive["call"].to_numpy()
full_oof_score = x_oof @ best_weights
summary = {
    "protocol": cfg["protocol"],
    "pool_members": pool_names,
    "pool_mean_outer_weights": {
        name: float(mean_weight[all_names.index(name)]) for name in pool_names
    },
    "selected_members": selected_names,
    "weights": {pool_names[i]: float(best_weights[i]) for i in selected_indices},
    "effective_member_count": float(1 / np.square(best_weights[selected_indices]).sum()),
    "honest_nested_metrics": {
        "AUROC": float(roc_auc_score(y, honest_score)),
        "AP": float(average_precision_score(y, honest_score)),
        "MCC_observed": float(matthews_corrcoef(y, honest_call)),
    },
    "descriptive_full_fit_metrics": {
        "AUROC": float(roc_auc_score(y, full_oof_score)),
        "AP": float(average_precision_score(y, full_oof_score)),
    },
    "test_positive_calls": int(candidate_test_call.sum()),
    "path": path,
    "candidate_created": False,
    "uploaded": False,
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(
    json.dumps(
        {key: value for key, value in summary.items() if key != "path"}, indent=2
    ),
    flush=True,
)
