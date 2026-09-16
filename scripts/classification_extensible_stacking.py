"""Nested cross-fit meta-models on an extensible classification prediction bank."""

import json
import os
from copy import deepcopy
from pathlib import Path

# These fits use many very small matrices. Large BLAS/PyTorch thread pools add
# orders of magnitude more scheduling overhead than useful work.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import polars as pl
import torch
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import rankdata
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, log_loss, matthews_corrcoef, roc_auc_score
from sklearn.preprocessing import StandardScaler

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
D6, A4 = "CYP2D6_is_TDI", "CYP3A4_is_TDI"
torch.set_num_threads(1)
torch.set_num_interop_threads(1)


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


def aligned_long(key: str, endpoint: str, method: str, ids: pl.DataFrame) -> np.ndarray:
    table = read(key)
    predicate = (pl.col("endpoint") == endpoint) & (pl.col("method") == method)
    if "strategy" in table.columns:
        predicate &= pl.col("strategy") == "umap"
    column = "probability" if "probability" in table.columns else "prediction"
    table = table.filter(predicate).select("compound_id", pl.col(column).alias("prediction"))
    joined = ids.join(table, on="compound_id", validate="1:1", maintain_order="left")
    if joined.height != ids.height or not joined["prediction"].is_finite().all():
        raise ValueError(f"Incomplete alignment for {key}/{endpoint}/{method}")
    return joined["prediction"].to_numpy()


def fixed_calls(values: np.ndarray, folds: np.ndarray, rate: float) -> np.ndarray:
    calls = np.zeros(len(values), dtype=bool)
    for fold in np.unique(folds):
        rows = np.flatnonzero(folds == fold)
        count = max(1, int(round(rate * len(rows))))
        calls[rows[np.argsort(-values[rows], kind="stable")[:count]]] = True
    return calls


def expected_mcc(y: np.ndarray, calls: np.ndarray, prevalence: float) -> float:
    positive = y == 1
    if not positive.any() or positive.all():
        return 0.0
    tpr, fpr = float(calls[positive].mean()), float(calls[~positive].mean())
    tp, fn = prevalence * tpr, prevalence * (1 - tpr)
    fp, tn = (1 - prevalence) * fpr, (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


def family(name: str) -> str:
    if name == "current_fixed_ratio":
        return "reference"
    if name.startswith("lgbm::"):
        return "tree"
    if name.startswith(("pfn::", "a4::tab")):
        return "tabular_foundation"
    if name.startswith("graph::"):
        return "graph"
    if name.startswith("morgan_") or name == "complement_nb":
        return "classical"
    if name.startswith("new::"):
        return next(item["family"] for item in cfg["candidate_families"] if item["name"] == name)
    return "other"


def build_bank(endpoint: str):
    if endpoint == D6:
        reference = (
            read("d6_reference_oof")
            .filter(
                (pl.col("endpoint") == endpoint) & (pl.col("method") == cfg["d6_reference_method"])
            )
            .sort("compound_id")
        )
        ids = reference.select("compound_id")
        folds = reference["fold"].to_numpy().astype(int)
        y = reference["y_true"].to_numpy().astype(int)
        groups = None
        specs = [
            ("lgbm::morgan_plain", "d6_lgbm_oof", "morgan_plain"),
            ("lgbm::mordred_plain", "d6_lgbm_oof", "mordred_plain"),
            ("lgbm::fp_mordred_plain", "d6_lgbm_oof", "fp_mordred_plain"),
            ("lgbm::external_class_plain", "d6_lgbm_oof", "external_class_plain"),
            ("lgbm::external_mixed_plain", "d6_lgbm_oof", "external_mixed_plain"),
            ("lgbm::external_class_mordred_plain", "d6_lgbm_oof", "external_class_mordred_plain"),
            ("pfn::external_class_pca256", "d6_pfn_oof", "external_class_pca256"),
            ("pfn::external_mixed_pca256", "d6_pfn_oof", "external_mixed_pca256"),
            ("pfn::chemeleon_pca256", "d6_pfn_oof", "chemeleon_pca256"),
            ("graph::tdi_only", "d6_graph_oof", "tdi_only"),
            ("graph::tdi_aux", "d6_graph_oof", "tdi_aux"),
        ]
        names = [name for name, _, _ in specs]
        values = [aligned_long(key, endpoint, method, ids) for name, key, method in specs]
        fixed = read("fixed_oof").filter(pl.col("endpoint") == endpoint)
        fixed_values = {}
        for name, method in [
            ("morgan_logistic", "logistic_c0.05"),
            ("morgan_rbf", "rbf_svm_c1"),
            ("complement_nb", "complement_nb_a1"),
        ]:
            table = fixed.filter(pl.col("method") == method).select(
                "compound_id", "prediction", "group"
            )
            joined = ids.join(table, on="compound_id", validate="1:1", maintain_order="left")
            fixed_values[name] = fold_rank(joined["prediction"].to_numpy(), folds)
            groups = joined["group"].cast(pl.String).to_numpy()
            names.append(name)
            values.append(joined["prediction"].to_numpy())
        old = fold_rank(reference["prediction"].to_numpy(), folds)
        live = (2 * old + fixed_values["morgan_logistic"]) / 3
        tails = {
            "tail_rbf": (2 * live + fixed_values["morgan_rbf"]) / 3,
            "tail_complement": (2 * live + fixed_values["complement_nb"]) / 3,
        }
        current = np.empty(len(y), dtype=float)
        for held, member in enumerate(cfg["d6_current_fold_members"]):
            current[folds == held] = tails[member][folds == held]
        raw = aligned_long("raw600_oof", endpoint, "external_mixed_raw600", ids)
        names.append("pfn::external_mixed_raw600")
        values.append(raw)
    else:
        table = read("a4_member_oof").sort("compound_id")
        ids = table.select("compound_id")
        folds = table["fold"].to_numpy().astype(int)
        y = table["y_true"].to_numpy().astype(int)
        base = [
            c for c in table.columns if c not in {"compound_id", "fold", "y_true", "caruana_auc"}
        ]
        names = [f"a4::{name}" for name in base]
        values = [table[name].to_numpy() for name in base]
        aux = aligned_long("a4_aux_oof", endpoint, cfg["a4_aux_method"], ids)
        names.append("a4::pxr_admetica")
        values.append(aux)
        current = 0.8 * fold_rank(table["caruana_auc"].to_numpy(), folds) + 0.2 * fold_rank(
            aux, folds
        )
        groups = None

    for item in cfg["candidate_families"]:
        method = item.get("method", item.get("method_by_endpoint", {}).get(endpoint))
        if method is not None:
            names.append(item["name"])
            values.append(aligned_long(item["source"], endpoint, method, ids))
            if groups is None and item["source"] == "mechanism_oof":
                source = (
                    read(item["source"])
                    .filter(
                        (pl.col("endpoint") == endpoint)
                        & (pl.col("strategy") == "umap")
                        & (pl.col("method") == method)
                    )
                    .select("compound_id", pl.col("group").cast(pl.String))
                )
                groups = ids.join(source, on="compound_id", validate="1:1", maintain_order="left")[
                    "group"
                ].to_numpy()
    names = ["current_fixed_ratio"] + names
    x = np.column_stack([fold_rank(current, folds)] + [fold_rank(value, folds) for value in values])
    if not np.isfinite(x).all() or groups is None:
        raise ValueError(f"Invalid bank for {endpoint}")
    return ids, folds, y, groups, names, x


def clipped_logloss(y: np.ndarray, p: np.ndarray) -> float:
    return float(log_loss(y, np.clip(p, 1e-6, 1 - 1e-6)))


def fit_ridge(x: np.ndarray, y: np.ndarray, c: float):
    scaler = StandardScaler().fit(x)
    model = LogisticRegression(C=c, max_iter=2000, solver="liblinear", random_state=20260915).fit(
        scaler.transform(x), y
    )
    return scaler, model


def ridge_predict(fitted, x: np.ndarray) -> np.ndarray:
    scaler, model = fitted
    return model.predict_proba(scaler.transform(x))[:, 1]


def simplex_loss(weights, x, y, loss_name, regularization):
    prediction = np.clip(x @ weights, 1e-6, 1 - 1e-6)
    if loss_name == "logloss":
        loss = clipped_logloss(y, prediction)
    elif loss_name == "brier":
        loss = float(np.mean((prediction - y) ** 2))
    else:
        raise ValueError("Pairwise loss uses its precomputed fit closure")
    return loss + regularization * float(np.sum(weights**2))


def fit_simplex(x, y, loss_name, regularization):
    count = x.shape[1]
    if loss_name == "pairwise":
        positive, negative = x[y == 1], x[y == 0]
        differences = (positive[:, None, :] - negative[None, :, :]).reshape(-1, count)
        temperature = cfg["pairwise_temperature"]

        def objective(weights):
            margins = differences @ weights / temperature
            ranking_loss = float(np.logaddexp(0, -margins).mean())
            return ranking_loss + regularization * float(np.sum(weights**2))

        def gradient(weights):
            margins = differences @ weights / temperature
            ranking_gradient = (-expit(-margins)[:, None] * differences / temperature).mean(axis=0)
            return ranking_gradient + 2 * regularization * weights

    else:

        def objective(weights):
            return simplex_loss(weights, x, y, loss_name, regularization)

        gradient = None

    result = minimize(
        objective,
        np.full(count, 1 / count),
        method="SLSQP",
        jac=gradient,
        bounds=[(0.0, 1.0)] * count,
        constraints={"type": "eq", "fun": lambda weights: weights.sum() - 1.0},
        options={"maxiter": 500, "ftol": 1e-10},
    )
    if not result.success:
        raise RuntimeError(result.message)
    weights = np.maximum(result.x, 0)
    return weights / weights.sum()


def choose_hyperparameter(x, y, folds, method, grid):
    scores = []
    for value in grid:
        pred = np.empty(len(y), dtype=float)
        for held in np.unique(folds):
            train, valid = folds != held, folds == held
            if method == "ridge":
                pred[valid] = ridge_predict(fit_ridge(x[train], y[train], value), x[valid])
            else:
                weights = fit_simplex(x[train], y[train], method, value)
                pred[valid] = x[valid] @ weights
        scores.append(clipped_logloss(y, pred))
    return grid[int(np.argmin(scores))], scores


class ExpertSet(torch.nn.Module):
    def __init__(self, token_width: int, hidden: int):
        super().__init__()
        self.phi = torch.nn.Sequential(
            torch.nn.Linear(token_width, hidden),
            torch.nn.Tanh(),
            torch.nn.Linear(hidden, hidden),
            torch.nn.Tanh(),
        )
        self.gate = torch.nn.Linear(hidden, 1)
        self.bias = torch.nn.Parameter(torch.tensor(0.0))
        self.raw_scale = torch.nn.Parameter(torch.tensor(1.0))

    def forward(self, tokens):
        gates = torch.softmax(self.gate(self.phi(tokens)).squeeze(-1), dim=1)
        pooled = torch.sum(gates * tokens[:, :, 0], dim=1).clamp(1e-5, 1 - 1e-5)
        scale = torch.nn.functional.softplus(self.raw_scale)
        logits = self.bias + scale * torch.logit(pooled)
        return logits, gates


def expert_tokens(
    x: np.ndarray, names: list[str], y: np.ndarray, train: np.ndarray, families: list[str]
):
    onehot = np.zeros((len(names), len(families)), dtype=np.float32)
    reliability = np.zeros((len(names), 3), dtype=np.float32)
    corr = np.corrcoef(x[train], rowvar=False)
    for index, name in enumerate(names):
        onehot[index, families.index(family(name))] = 1
        reliability[index] = [
            roc_auc_score(y[train], x[train, index]),
            average_precision_score(y[train], x[train, index]),
            1 - (np.abs(corr[index]).sum() - 1) / max(1, len(names) - 1),
        ]
    static = np.column_stack([reliability, onehot])
    static[:, :3] = (static[:, :3] - static[:, :3].mean(axis=0)) / (
        static[:, :3].std(axis=0) + 1e-6
    )
    return np.concatenate(
        [x[:, :, None].astype(np.float32), np.broadcast_to(static, (len(x), *static.shape))],
        axis=2,
    )


def fit_expert_set(tokens, y, train, valid, seed):
    settings = cfg["deepsets"]
    torch.manual_seed(seed)
    model = ExpertSet(tokens.shape[2], settings["hidden"])
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=settings["learning_rate"], weight_decay=settings["weight_decay"]
    )
    xt = torch.from_numpy(tokens)
    yt = torch.from_numpy(y.astype(np.float32))
    best_loss, best_state, patience = np.inf, None, 0
    for _ in range(settings["max_epochs"]):
        model.train()
        optimizer.zero_grad()
        logits, _ = model(xt[train])
        loss = torch.nn.functional.binary_cross_entropy_with_logits(logits, yt[train])
        loss.backward()
        optimizer.step()
        model.eval()
        with torch.no_grad():
            valid_logits, _ = model(xt[valid])
            value = float(
                torch.nn.functional.binary_cross_entropy_with_logits(valid_logits, yt[valid])
            )
        if value < best_loss - 1e-5:
            best_loss, best_state, patience = value, deepcopy(model.state_dict()), 0
        else:
            patience += 1
            if patience >= settings["patience"]:
                break
    model.load_state_dict(best_state)
    return model


def crossfit(endpoint, folds, y, names, x):
    methods = ["ridge", "simplex_logloss", "simplex_brier", "simplex_pairwise", "expert_set"]
    predictions = {name: np.empty(len(y), dtype=float) for name in methods}
    fits = []
    all_families = sorted({family(name) for name in names})
    for held in np.unique(folds):
        outer_train, outer_valid = folds != held, folds == held
        local_folds = folds[outer_train]
        c, c_scores = choose_hyperparameter(
            x[outer_train], y[outer_train], local_folds, "ridge", cfg["ridge_c"]
        )
        predictions["ridge"][outer_valid] = ridge_predict(
            fit_ridge(x[outer_train], y[outer_train], c), x[outer_valid]
        )
        fits.append(
            {
                "endpoint": endpoint,
                "method": "ridge",
                "held_fold": int(held),
                "hyperparameter": c,
                "inner_logloss": json.dumps(c_scores),
                "weights": None,
            }
        )
        for loss_name in ["logloss", "brier", "pairwise"]:
            value, scores = choose_hyperparameter(
                x[outer_train], y[outer_train], local_folds, loss_name, cfg["simplex_lambda"]
            )
            weights = fit_simplex(x[outer_train], y[outer_train], loss_name, value)
            key = f"simplex_{loss_name}"
            predictions[key][outer_valid] = x[outer_valid] @ weights
            fits.append(
                {
                    "endpoint": endpoint,
                    "method": key,
                    "held_fold": int(held),
                    "hyperparameter": value,
                    "inner_logloss": json.dumps(scores),
                    "weights": json.dumps(dict(zip(names, weights.tolist(), strict=True))),
                }
            )
        deep_predictions = []
        gate_values = []
        outer_rows = np.flatnonzero(outer_train)
        for inner_held in np.unique(local_folds):
            inner_train = outer_rows[local_folds != inner_held]
            inner_valid = outer_rows[local_folds == inner_held]
            tokens = expert_tokens(x, names, y, inner_train, all_families)
            for seed in cfg["deepsets"]["seeds"]:
                model = fit_expert_set(tokens, y, inner_train, inner_valid, seed + int(held) * 100)
                model.eval()
                with torch.no_grad():
                    logits, gates = model(torch.from_numpy(tokens[outer_valid]))
                deep_predictions.append(expit(logits.numpy()))
                gate_values.append(gates.numpy().mean(axis=0))
        predictions["expert_set"][outer_valid] = np.mean(deep_predictions, axis=0)
        mean_gates = np.mean(gate_values, axis=0)
        fits.append(
            {
                "endpoint": endpoint,
                "method": "expert_set",
                "held_fold": int(held),
                "hyperparameter": None,
                "inner_logloss": None,
                "weights": json.dumps(dict(zip(names, mean_gates.tolist(), strict=True))),
            }
        )
    return predictions, fits


metric_rows, fit_rows, prediction_frames, bootstrap_rows, bank_frames = [], [], [], [], []
rng = np.random.default_rng(cfg["bootstrap_seed"])
for endpoint in [D6, A4]:
    ids, folds, y, groups, names, x = build_bank(endpoint)
    settings = cfg["endpoint_settings"][endpoint]
    predictions = {"current_fixed_ratio": x[:, 0], "equal_mean_all": x.mean(axis=1)}
    stacked, fits = crossfit(endpoint, folds, y, names, x)
    predictions.update(stacked)
    fit_rows.extend(fits)
    for index, name in enumerate(names):
        bank_frames.append(
            ids.with_columns(
                pl.lit(endpoint).alias("endpoint"),
                pl.Series("fold", folds),
                pl.Series("y_true", y),
                pl.lit(name).alias("model"),
                pl.Series("prediction", x[:, index]),
            )
        )
    calls_by_method = {}
    for method, prediction in predictions.items():
        calls = fixed_calls(prediction, folds, settings["call_rate"])
        calls_by_method[method] = calls
        metric_rows.append(
            {
                "endpoint": endpoint,
                "method": method,
                "experts": len(names),
                "AUROC": float(roc_auc_score(y, prediction)),
                "AP": float(average_precision_score(y, prediction)),
                "logloss": clipped_logloss(y, prediction),
                "MCC_observed": float(matthews_corrcoef(y, calls)),
                "MCC_target_prior": expected_mcc(y, calls, settings["target_prior"]),
            }
        )
        prediction_frames.append(
            ids.with_columns(
                pl.lit(endpoint).alias("endpoint"),
                pl.Series("fold", folds),
                pl.Series("y_true", y),
                pl.lit(method).alias("method"),
                pl.Series("prediction", prediction),
                pl.Series("call", calls),
            )
        )
    group_values = np.unique(groups)
    group_rows = {group: np.flatnonzero(groups == group) for group in group_values}
    for method in predictions:
        if method == "current_fixed_ratio":
            continue
        deltas = []
        for _ in range(cfg["bootstrap_replicates"]):
            sampled = rng.choice(group_values, size=len(group_values), replace=True)
            rows = np.concatenate([group_rows[group] for group in sampled])
            deltas.append(
                expected_mcc(y[rows], calls_by_method[method][rows], settings["target_prior"])
                - expected_mcc(
                    y[rows], calls_by_method["current_fixed_ratio"][rows], settings["target_prior"]
                )
            )
        delta = np.asarray(deltas)
        bootstrap_rows.append(
            {
                "endpoint": endpoint,
                "method": method,
                "mean_MCC_delta": float(delta.mean()),
                "ci025": float(np.quantile(delta, 0.025)),
                "ci975": float(np.quantile(delta, 0.975)),
                "probability_positive": float((delta > 0).mean()),
            }
        )

metrics = pl.DataFrame(metric_rows).sort(["endpoint", "MCC_target_prior"], descending=[False, True])
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(bootstrap_rows).write_csv(out / "bootstrap.csv")
pl.DataFrame(fit_rows).write_parquet(out / "fits.parquet")
pl.concat(prediction_frames).write_parquet(out / "oof.parquet")
pl.concat(bank_frames).write_parquet(out / "prediction-bank.parquet")
(out / "summary.json").write_text(json.dumps({"protocol": cfg["protocol"]}, indent=2) + "\n")
print(metrics, flush=True)
