"""Cross-fit live classification pools with representation peers and replacements."""

import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import rankdata, spearmanr
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
D6, A4 = "CYP2D6_is_TDI", "CYP3A4_is_TDI"


def read(name: str) -> pl.DataFrame:
    return pl.read_parquet(inp / cfg[name])


def unit_rank(values: np.ndarray) -> np.ndarray:
    return rankdata(values, method="average") / (len(values) + 1.0)


def model_score(model: object, features: np.ndarray) -> np.ndarray:
    if hasattr(model, "predict_proba"):
        return model.predict_proba(features)[:, 1]
    return model.decision_function(features)


def fold_rank(values: np.ndarray, folds: np.ndarray) -> np.ndarray:
    result = np.empty(len(values), dtype=float)
    for fold in np.unique(folds):
        mask = folds == fold
        result[mask] = unit_rank(values[mask])
    return result


def fixed_calls(values: np.ndarray, rate: float, folds: np.ndarray) -> np.ndarray:
    result = np.zeros(len(values), dtype=bool)
    for fold in np.unique(folds):
        rows = np.flatnonzero(folds == fold)
        count = max(1, int(round(rate * len(rows))))
        result[rows[np.argsort(-values[rows], kind="stable")[:count]]] = True
    return result


def top_calls(values: np.ndarray, count: int) -> np.ndarray:
    result = np.zeros(len(values), dtype=bool)
    result[np.argsort(-values, kind="stable")[:count]] = True
    return result


def expected_mcc(y: np.ndarray, selected: np.ndarray, prevalence: float) -> float:
    positive = y == 1
    if not positive.any() or positive.all():
        return 0.0
    tpr, fpr = float(selected[positive].mean()), float(selected[~positive].mean())
    tp, fn = prevalence * tpr, prevalence * (1 - tpr)
    fp, tn = (1 - prevalence) * fpr, (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


def objective(
    y: np.ndarray,
    prediction: np.ndarray,
    folds: np.ndarray,
    recipe: str,
    rate: float,
    prevalence: float,
) -> float:
    if recipe == "caruana_auc":
        return float(roc_auc_score(y, prediction))
    return expected_mcc(y, fixed_calls(prediction, rate, folds), prevalence)


def fit_caruana(
    x: np.ndarray,
    y: np.ndarray,
    folds: np.ndarray,
    names: list[str],
    recipe: str,
    rate: float,
    prevalence: float,
) -> tuple[np.ndarray, list[dict]]:
    cap = 5 if recipe.endswith("sparse5") else len(names)
    counts = np.zeros(len(names), dtype=int)
    total = np.zeros(len(y), dtype=float)
    best_value, best_counts = -np.inf, None
    history = []
    for step in range(1, cfg["steps"] + 1):
        candidates = [
            index for index in range(len(names)) if counts[index] or np.count_nonzero(counts) < cap
        ]
        values = [
            objective(
                y,
                (total + x[:, index]) / step,
                folds,
                recipe,
                rate,
                prevalence,
            )
            for index in candidates
        ]
        choice = candidates[int(np.argmax(values))]
        total += x[:, choice]
        counts[choice] += 1
        value = max(values)
        history.append({"step": step, "value": value, "added": names[choice]})
        if value > best_value + 1e-12:
            best_value, best_counts = value, counts.copy()
    if best_counts is None:
        raise RuntimeError("Caruana selection did not produce weights")
    return best_counts / best_counts.sum(), history


def aligned_long(
    source: str,
    endpoint: str,
    method: str,
    reference_ids: pl.DataFrame,
    test: bool = False,
) -> np.ndarray:
    table = read(source)
    predicate = (pl.col("endpoint") == endpoint) & (pl.col("method") == method)
    if "strategy" in table.columns:
        predicate &= pl.col("strategy") == "umap"
    table = table.filter(predicate)
    if test and table.select("compound_id").is_duplicated().any():
        table = table.group_by("compound_id").agg(pl.col("prediction").mean())
    table = table.select("compound_id", "prediction")
    joined = reference_ids.join(table, on="compound_id", validate="1:1", maintain_order="left")
    if joined.height != reference_ids.height or not joined["prediction"].is_finite().all():
        raise ValueError(f"Incomplete alignment for {source}/{endpoint}/{method}")
    return joined["prediction"].to_numpy()


def aligned_d6_atomic(
    source: str,
    method: str,
    prediction_column: str,
    reference_ids: pl.DataFrame,
    test: bool = False,
) -> np.ndarray:
    table = read(source).filter((pl.col("endpoint") == D6) & (pl.col("method") == method))
    if test:
        table = table.group_by("compound_id").agg(
            pl.col(prediction_column).mean().alias("prediction")
        )
    else:
        table = table.select("compound_id", pl.col(prediction_column).alias("prediction"))
    joined = reference_ids.join(table, on="compound_id", validate="1:1", maintain_order="left")
    if joined.height != reference_ids.height or not joined["prediction"].is_finite().all():
        raise ValueError(f"Incomplete D6 atomic alignment for {source}/{method}")
    return joined["prediction"].to_numpy()


def metric_row(
    endpoint: str,
    pool: str,
    recipe: str,
    y: np.ndarray,
    prediction: np.ndarray,
    calls: np.ndarray,
    prevalence: float,
    test_prediction: np.ndarray,
    current_test: np.ndarray,
    test_count: int,
) -> dict:
    return {
        "endpoint": endpoint,
        "pool": pool,
        "recipe": recipe,
        "AUROC": float(roc_auc_score(y, prediction)),
        "AP": float(average_precision_score(y, prediction)),
        "MCC_observed": float(matthews_corrcoef(y, calls)),
        "MCC_target_prior": expected_mcc(y, calls, prevalence),
        "calls": int(calls.sum()),
        "test_spearman_to_current": float(spearmanr(test_prediction, current_test).statistic),
        "test_call_swaps_each_way": int(
            (top_calls(test_prediction, test_count) & ~top_calls(current_test, test_count)).sum()
        ),
    }


def crossfit_pool(
    endpoint: str,
    pool_name: str,
    recipe: str,
    names: list[str],
    x: np.ndarray,
    test_x: np.ndarray,
    y: np.ndarray,
    folds: np.ndarray,
    rate: float,
    prevalence: float,
    selection_rows: list[dict],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    prediction = np.empty(len(y), dtype=float)
    for held in np.unique(folds):
        train, valid = folds != held, folds == held
        weights, history = fit_caruana(
            x[train], y[train], folds[train], names, recipe, rate, prevalence
        )
        prediction[valid] = x[valid] @ weights
        selection_rows.append(
            {
                "endpoint": endpoint,
                "pool": pool_name,
                "recipe": recipe,
                "held_fold": int(held),
                "weights": json.dumps(dict(zip(names, weights.tolist(), strict=True))),
                "history": json.dumps(history),
            }
        )
    weights, history = fit_caruana(x, y, folds, names, recipe, rate, prevalence)
    selection_rows.append(
        {
            "endpoint": endpoint,
            "pool": pool_name,
            "recipe": recipe,
            "held_fold": -1,
            "weights": json.dumps(dict(zip(names, weights.tolist(), strict=True))),
            "history": json.dumps(history),
        }
    )
    return prediction, test_x @ weights, weights


# Common A4 identity, labels, groups, live top-level members and exact current score.
a4_oof = read("a4_member_oof").sort("compound_id")
a4_test = read("a4_member_test").sort("compound_id")
a4_ids, a4_test_ids = a4_oof.select("compound_id"), a4_test.select("compound_id")
a4_folds = a4_oof["fold"].to_numpy().astype(int)
a4_y = a4_oof["y_true"].to_numpy().astype(int)
a4_base_names = [
    column
    for column in a4_oof.columns
    if column not in {"compound_id", "fold", "y_true", "caruana_auc"}
]
a4_base_oof = [fold_rank(a4_oof[name].to_numpy(), a4_folds) for name in a4_base_names]
a4_base_test = [unit_rank(a4_test[name].to_numpy()) for name in a4_base_names]
a4_aux_oof = aligned_long("a4_aux_oof", A4, cfg["a4_aux_method"], a4_ids)
a4_aux_test = aligned_long("a4_aux_test", A4, cfg["a4_aux_method"], a4_test_ids, True)
a4_base_names.append("pxr_admetica")
a4_base_oof.append(fold_rank(a4_aux_oof, a4_folds))
a4_base_test.append(unit_rank(a4_aux_test))
a4_current = 0.8 * fold_rank(a4_oof["caruana_auc"].to_numpy(), a4_folds) + 0.2 * fold_rank(
    a4_aux_oof, a4_folds
)
a4_current_test = 0.8 * unit_rank(a4_test["caruana_auc"].to_numpy()) + 0.2 * unit_rank(a4_aux_test)

# D6 live top-level members and exact current score.
d6_oof = (
    read("d6_reference_oof")
    .filter((pl.col("endpoint") == D6) & (pl.col("method") == cfg["d6_reference_method"]))
    .sort("compound_id")
)
d6_test = (
    read("d6_reference_test")
    .filter((pl.col("endpoint") == D6) & (pl.col("method") == cfg["d6_reference_method"]))
    .sort("compound_id")
)
d6_ids, d6_test_ids = d6_oof.select("compound_id"), d6_test.select("compound_id")
d6_folds = d6_oof["fold"].to_numpy().astype(int)
d6_y = d6_oof["y_true"].to_numpy().astype(int)
d6_old = fold_rank(d6_oof["prediction"].to_numpy(), d6_folds)
d6_old_test = unit_rank(d6_test["prediction"].to_numpy())
d6_fixed = read("fixed_oof").sort("compound_id")
d6_fixed_oof, d6_fixed_test = {}, {}

compounds = read("compounds").sort("compound_id")
test_mapping = read("test_mapping").select("compound_id").unique().sort("compound_id")
if not test_mapping["compound_id"].equals(d6_test["compound_id"]):
    raise ValueError("D6 test identities differ")
id_to_index = {compound_id: index for index, compound_id in enumerate(compounds["compound_id"])}
test_indices = np.asarray([id_to_index[compound_id] for compound_id in test_mapping["compound_id"]])
generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
morgan = np.stack(
    [
        generator.GetFingerprintAsNumPy(Chem.MolFromSmiles(smiles))
        for smiles in compounds["standardized_smiles"]
    ]
).astype(np.float32)
for method, model_key in [
    ("logistic_c0.05", "logistic_models"),
    ("complement_nb_a1", "complement_models"),
    ("rbf_svm_c1", "rbf_models"),
]:
    table = d6_fixed.filter((pl.col("endpoint") == D6) & (pl.col("method") == method)).sort(
        "compound_id"
    )
    if not table["compound_id"].equals(d6_oof["compound_id"]):
        raise ValueError(f"D6 fixed OOF identities differ for {method}")
    d6_fixed_oof[method] = fold_rank(table["prediction"].to_numpy(), d6_folds)
    fold_test = []
    for model_file in cfg[model_key]:
        model = joblib.load(inp / model_file)
        fold_test.append(unit_rank(model_score(model, morgan[test_indices])))
    d6_fixed_test[method] = unit_rank(np.mean(fold_test, axis=0))
d6_live = (2 * d6_old + d6_fixed_oof["logistic_c0.05"]) / 3
d6_live_test = (2 * d6_old_test + d6_fixed_test["logistic_c0.05"]) / 3
d6_variants_oof = {
    "live": d6_live,
    "live2_rbf": (2 * d6_live + d6_fixed_oof["rbf_svm_c1"]) / 3,
    "live2_complement": (2 * d6_live + d6_fixed_oof["complement_nb_a1"]) / 3,
}
d6_variants_test = {
    "live": d6_live_test,
    "live2_rbf": (2 * d6_live_test + d6_fixed_test["rbf_svm_c1"]) / 3,
    "live2_complement": (2 * d6_live_test + d6_fixed_test["complement_nb_a1"]) / 3,
}
d6_current = np.empty(len(d6_y), dtype=float)
for fold, member in enumerate(cfg["d6_current_fold_members"]):
    mask = d6_folds == fold
    d6_current[mask] = d6_variants_oof[member][mask]
d6_current_test = d6_variants_test[cfg["d6_deployment_member"]]
d6_base_names = list(d6_variants_oof)
d6_base_oof = [d6_variants_oof[name] for name in d6_base_names]
d6_base_test = [d6_variants_test[name] for name in d6_base_names]

# A genuinely flat legacy pool: no old11 mean and no live/live2 wrappers.
d6_atomic_specs = [
    (
        "lgbm::morgan_plain",
        "d6_lgbm_oof",
        "d6_lgbm_test",
        "morgan_plain",
        "prediction",
    ),
    (
        "lgbm::mordred_plain",
        "d6_lgbm_oof",
        "d6_lgbm_test",
        "mordred_plain",
        "prediction",
    ),
    (
        "lgbm::fp_mordred_plain",
        "d6_lgbm_oof",
        "d6_lgbm_test",
        "fp_mordred_plain",
        "prediction",
    ),
    (
        "lgbm::external_class_plain",
        "d6_lgbm_oof",
        "d6_lgbm_test",
        "external_class_plain",
        "prediction",
    ),
    (
        "lgbm::external_mixed_plain",
        "d6_lgbm_oof",
        "d6_lgbm_test",
        "external_mixed_plain",
        "prediction",
    ),
    (
        "lgbm::external_class_mordred_plain",
        "d6_lgbm_oof",
        "d6_lgbm_test",
        "external_class_mordred_plain",
        "prediction",
    ),
    (
        "pfn::external_class_pca256",
        "d6_pfn_oof",
        "d6_pfn_test",
        "external_class_pca256",
        "probability",
    ),
    (
        "pfn::external_mixed_pca256",
        "d6_pfn_oof",
        "d6_pfn_test",
        "external_mixed_pca256",
        "probability",
    ),
    (
        "pfn::chemeleon_pca256",
        "d6_pfn_oof",
        "d6_pfn_test",
        "chemeleon_pca256",
        "probability",
    ),
    (
        "graph::tdi_only",
        "d6_graph_oof",
        "d6_graph_test",
        "tdi_only",
        "probability",
    ),
    (
        "graph::tdi_aux",
        "d6_graph_oof",
        "d6_graph_test",
        "tdi_aux",
        "probability",
    ),
]
d6_flat_names = [spec[0] for spec in d6_atomic_specs] + [
    "morgan_logistic",
    "morgan_rbf",
    "complement_nb",
]
d6_flat_oof = [
    fold_rank(aligned_d6_atomic(oof, method, column, d6_ids), d6_folds)
    for _, oof, _, method, column in d6_atomic_specs
] + [
    d6_fixed_oof["logistic_c0.05"],
    d6_fixed_oof["rbf_svm_c1"],
    d6_fixed_oof["complement_nb_a1"],
]
d6_flat_test = [
    unit_rank(aligned_d6_atomic(test, method, column, d6_test_ids, True))
    for _, _, test, method, column in d6_atomic_specs
] + [
    d6_fixed_test["logistic_c0.05"],
    d6_fixed_test["rbf_svm_c1"],
    d6_fixed_test["complement_nb_a1"],
]

# Load one prespecified representative for each previously screened family.
source_keys = {
    "mechanism": ("mechanism_oof", "mechanism_test"),
    "activity": ("activity_oof", "activity_test"),
    "mace": ("mace_oof", "mace_test"),
    "embedding": ("embedding_oof", "embedding_test"),
}
endpoint_data = {
    D6: {
        "ids": d6_ids,
        "test_ids": d6_test_ids,
        "folds": d6_folds,
        "y": d6_y,
        "base_names": d6_base_names,
        "base_oof": d6_base_oof,
        "base_test": d6_base_test,
        "current": d6_current,
        "current_test": d6_current_test,
    },
    A4: {
        "ids": a4_ids,
        "test_ids": a4_test_ids,
        "folds": a4_folds,
        "y": a4_y,
        "base_names": a4_base_names,
        "base_oof": a4_base_oof,
        "base_test": a4_base_test,
        "current": a4_current,
        "current_test": a4_current_test,
    },
}
for endpoint, data in endpoint_data.items():
    new_names, new_oof, new_test = [], [], []
    raw_new_oof, raw_new_test = [], []
    for family in cfg["candidate_families"]:
        method = family.get("method", family.get("method_by_endpoint", {}).get(endpoint))
        if method is None:
            continue
        oof_key, test_key = source_keys[family["source"]]
        new_names.append(family["name"])
        candidate_oof = fold_rank(
            aligned_long(oof_key, endpoint, method, data["ids"]), data["folds"]
        )
        candidate_test = unit_rank(aligned_long(test_key, endpoint, method, data["test_ids"], True))
        raw_new_oof.append(candidate_oof)
        raw_new_test.append(candidate_test)
        if endpoint == D6:
            candidate_oof = (2 * d6_live + candidate_oof) / 3
            candidate_test = (2 * d6_live_test + candidate_test) / 3
        new_oof.append(candidate_oof)
        new_test.append(candidate_test)
    data["new_names"], data["new_oof"], data["new_test"] = new_names, new_oof, new_test
    data["raw_new_oof"], data["raw_new_test"] = raw_new_oof, raw_new_test

# Cross-fit base, expanded, and explicit base-minus-one-plus-one replacement pools.
metric_rows, selection_rows, prediction_rows, test_rows = [], [], [], []
bootstrap_rows = []
rng = np.random.default_rng(cfg["bootstrap_seed"])
for endpoint, data in endpoint_data.items():
    settings = cfg["endpoint_settings"][endpoint]
    folds, y = data["folds"], data["y"]
    group_source = (
        read("mechanism_oof")
        .filter(
            (pl.col("endpoint") == endpoint)
            & (pl.col("strategy") == "umap")
            & (pl.col("method") == cfg["candidate_families"][0]["method"])
        )
        .select("compound_id", "group")
    )
    groups = (
        data["ids"]
        .join(
            group_source,
            on="compound_id",
            validate="1:1",
            maintain_order="left",
        )["group"]
        .cast(pl.String)
        .to_numpy()
    )
    base_x, base_test_x = np.column_stack(data["base_oof"]), np.column_stack(data["base_test"])
    new_x, new_test_x = np.column_stack(data["new_oof"]), np.column_stack(data["new_test"])
    pools = {
        "base": (data["base_names"], base_x, base_test_x),
        "expanded": (
            data["base_names"] + data["new_names"],
            np.column_stack([base_x, new_x]),
            np.column_stack([base_test_x, new_test_x]),
        ),
    }
    if endpoint == A4:
        for pool_name, member_names in cfg["a4_compressed_pools"].items():
            member_indices = [data["base_names"].index(name) for name in member_names]
            pools[pool_name] = (
                member_names,
                base_x[:, member_indices],
                base_test_x[:, member_indices],
            )
    else:
        flat_x = np.column_stack(d6_flat_oof)
        flat_test_x = np.column_stack(d6_flat_test)
        raw_new_x = np.column_stack(data["raw_new_oof"])
        raw_new_test_x = np.column_stack(data["raw_new_test"])
        pools["flat_legacy"] = (d6_flat_names, flat_x, flat_test_x)
        pools["flat_expanded"] = (
            d6_flat_names + data["new_names"],
            np.column_stack([flat_x, raw_new_x]),
            np.column_stack([flat_test_x, raw_new_test_x]),
        )
    current_calls = fixed_calls(data["current"], settings["call_rate"], folds)
    metric_rows.append(
        metric_row(
            endpoint,
            "current_fixed",
            "current_fixed",
            y,
            data["current"],
            current_calls,
            settings["target_prior"],
            data["current_test"],
            data["current_test"],
            settings["test_calls"],
        )
    )
    for recipe in cfg["recipes"]:
        evaluated = {}
        for pool_name, (names, x, test_x) in pools.items():
            pred, test_pred, _ = crossfit_pool(
                endpoint,
                pool_name,
                recipe,
                names,
                x,
                test_x,
                y,
                folds,
                settings["call_rate"],
                settings["target_prior"],
                selection_rows,
            )
            calls = fixed_calls(pred, settings["call_rate"], folds)
            metric_rows.append(
                metric_row(
                    endpoint,
                    pool_name,
                    recipe,
                    y,
                    pred,
                    calls,
                    settings["target_prior"],
                    test_pred,
                    data["current_test"],
                    settings["test_calls"],
                )
            )
            evaluated[pool_name] = (pred, calls)
            for index, compound_id in enumerate(data["ids"]["compound_id"]):
                prediction_rows.append(
                    {
                        "compound_id": compound_id,
                        "endpoint": endpoint,
                        "fold": int(folds[index]),
                        "group": groups[index],
                        "y_true": int(y[index]),
                        "pool": pool_name,
                        "recipe": recipe,
                        "prediction": float(pred[index]),
                        "call": bool(calls[index]),
                    }
                )
            for index, compound_id in enumerate(data["test_ids"]["compound_id"]):
                test_rows.append(
                    {
                        "compound_id": compound_id,
                        "endpoint": endpoint,
                        "pool": pool_name,
                        "recipe": recipe,
                        "prediction": float(test_pred[index]),
                        "call": bool(top_calls(test_pred, settings["test_calls"])[index]),
                    }
                )
        group_values = np.unique(groups)
        group_rows = {group: np.flatnonzero(groups == group) for group in group_values}
        base_pred, base_calls = evaluated["base"]
        for contrast_pool, (contrast_pred, contrast_calls) in evaluated.items():
            if contrast_pool == "base":
                continue
            reference_name = "flat_legacy" if contrast_pool == "flat_expanded" else "base"
            reference_pred, reference_calls = evaluated[reference_name]
            delta_mcc, delta_auc = [], []
            for _ in range(cfg["bootstrap_replicates"]):
                sampled = rng.choice(group_values, size=len(group_values), replace=True)
                rows = np.concatenate([group_rows[group] for group in sampled])
                delta_mcc.append(
                    expected_mcc(y[rows], contrast_calls[rows], settings["target_prior"])
                    - expected_mcc(y[rows], reference_calls[rows], settings["target_prior"])
                )
                if np.unique(y[rows]).size == 2:
                    delta_auc.append(
                        roc_auc_score(y[rows], contrast_pred[rows])
                        - roc_auc_score(y[rows], reference_pred[rows])
                    )
            for metric, values in [
                ("MCC_target_prior", delta_mcc),
                ("AUROC", delta_auc),
            ]:
                values = np.asarray(values)
                bootstrap_rows.append(
                    {
                        "endpoint": endpoint,
                        "recipe": recipe,
                        "contrast": f"{contrast_pool}-minus-{reference_name}",
                        "metric": metric,
                        "mean": float(values.mean()),
                        "ci025": float(np.quantile(values, 0.025)),
                        "ci975": float(np.quantile(values, 0.975)),
                        "probability_positive": float((values > 0).mean()),
                    }
                )

        # Explicitly reselect after replacing one incumbent by one new family.
        for removed_index, removed in enumerate(data["base_names"]):
            keep = [index for index in range(len(data["base_names"])) if index != removed_index]
            for new_index, added in enumerate(data["new_names"]):
                names = [data["base_names"][index] for index in keep] + [added]
                x = np.column_stack([base_x[:, keep], new_x[:, new_index]])
                test_x = np.column_stack([base_test_x[:, keep], new_test_x[:, new_index]])
                pool_name = f"swap::{removed}=>{added}"
                pred, test_pred, _ = crossfit_pool(
                    endpoint,
                    pool_name,
                    recipe,
                    names,
                    x,
                    test_x,
                    y,
                    folds,
                    settings["call_rate"],
                    settings["target_prior"],
                    selection_rows,
                )
                calls = fixed_calls(pred, settings["call_rate"], folds)
                row = metric_row(
                    endpoint,
                    pool_name,
                    recipe,
                    y,
                    pred,
                    calls,
                    settings["target_prior"],
                    test_pred,
                    data["current_test"],
                    settings["test_calls"],
                )
                row.update({"removed": removed, "added": added})
                metric_rows.append(row)
                for index, compound_id in enumerate(data["ids"]["compound_id"]):
                    prediction_rows.append(
                        {
                            "compound_id": compound_id,
                            "endpoint": endpoint,
                            "fold": int(folds[index]),
                            "group": groups[index],
                            "y_true": int(y[index]),
                            "pool": pool_name,
                            "recipe": recipe,
                            "prediction": float(pred[index]),
                            "call": bool(calls[index]),
                        }
                    )
                test_calls = top_calls(test_pred, settings["test_calls"])
                for index, compound_id in enumerate(data["test_ids"]["compound_id"]):
                    test_rows.append(
                        {
                            "compound_id": compound_id,
                            "endpoint": endpoint,
                            "pool": pool_name,
                            "recipe": recipe,
                            "prediction": float(test_pred[index]),
                            "call": bool(test_calls[index]),
                        }
                    )
                deltas = []
                for _ in range(cfg["bootstrap_replicates"]):
                    sampled = rng.choice(group_values, size=len(group_values), replace=True)
                    rows = np.concatenate([group_rows[group] for group in sampled])
                    deltas.append(
                        expected_mcc(y[rows], calls[rows], settings["target_prior"])
                        - expected_mcc(y[rows], base_calls[rows], settings["target_prior"])
                    )
                values = np.asarray(deltas)
                bootstrap_rows.append(
                    {
                        "endpoint": endpoint,
                        "recipe": recipe,
                        "contrast": f"{pool_name}-minus-base",
                        "metric": "MCC_target_prior",
                        "mean": float(values.mean()),
                        "ci025": float(np.quantile(values, 0.025)),
                        "ci975": float(np.quantile(values, 0.975)),
                        "probability_positive": float((values > 0).mean()),
                    }
                )

pl.DataFrame(metric_rows).write_csv(out / "metrics.csv")
pl.DataFrame(selection_rows).write_parquet(out / "selections.parquet")
pl.DataFrame(prediction_rows).write_parquet(out / "oof.parquet")
pl.DataFrame(test_rows).write_parquet(out / "test.parquet")
pl.DataFrame(bootstrap_rows).write_csv(out / "bootstrap.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "selection": cfg["selection"],
            "candidate_families": cfg["candidate_families"],
            "caveats": [
                (
                    "This is second-level cross-fit over previously generated OOF predictions, "
                    "not fresh nested base-model refitting."
                ),
                (
                    "Candidate families were screened earlier; expanded-pool results remain "
                    "conditional on that history."
                ),
                (
                    "One-for-one swap screens are diagnostic multiplicity, not independent "
                    "model selection evidence."
                ),
            ],
        },
        indent=2,
    )
    + "\n"
)
print(
    pl.DataFrame(metric_rows).filter(pl.col("pool").is_in(["current_fixed", "base", "expanded"])),
    flush=True,
)
