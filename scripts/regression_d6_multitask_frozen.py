"""Four-head neural readout with weak-label transfer and a D6-specific head."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from torch import nn
from torch.nn import functional as F


def st_rae(frame: pl.DataFrame) -> float:
    y = frame["y_true"].to_numpy()
    pred = frame["prediction"].to_numpy()
    lo = frame["conf_low"].to_numpy()
    hi = frame["conf_high"].to_numpy()
    center = y.mean()
    return float(
        (np.maximum(lo - pred, 0) + np.maximum(pred - hi, 0)).sum()
        / (np.maximum(lo - center, 0) + np.maximum(center - hi, 0)).sum()
    )


def make_model(input_dim: int, heads: int, cfg: dict) -> nn.Module:
    h1, h2 = cfg["hidden"]
    return nn.Sequential(
        nn.Linear(input_dim, h1),
        nn.GELU(),
        nn.Dropout(cfg["dropout"]),
        nn.Linear(h1, h2),
        nn.GELU(),
        nn.Dropout(cfg["dropout"]),
        nn.Linear(h2, heads),
    )


def fit_predict(
    x: np.ndarray,
    y: np.ndarray,
    primary: np.ndarray,
    train_idx: np.ndarray,
    held_idx: np.ndarray,
    arm: dict,
    cfg: dict,
    weak_cutoff: float,
    d6_col: int,
    seed: int,
    device: torch.device,
) -> tuple[np.ndarray, dict]:
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    train_y = y[train_idx].copy()
    mask = np.isfinite(train_y)
    if arm["donor_labels"] == "none":
        mask[:, np.arange(y.shape[1]) != d6_col] = False
    elif arm["donor_labels"] == "weak":
        donor_cols = np.arange(y.shape[1]) != d6_col
        mask[:, donor_cols] &= primary[train_idx][:, donor_cols] >= weak_cutoff
    elif arm["donor_labels"] != "all":
        raise ValueError(arm)
    if mask[:, d6_col].sum() < 1000:
        raise ValueError("Insufficient D6 training labels")
    means = np.zeros(y.shape[1], dtype=np.float32)
    stds = np.ones(y.shape[1], dtype=np.float32)
    for col in range(y.shape[1]):
        if mask[:, col].any():
            means[col] = train_y[mask[:, col], col].mean()
            stds[col] = max(train_y[mask[:, col], col].std(), 0.2)
    normalized = np.where(mask, (np.nan_to_num(train_y) - means) / stds, 0).astype(np.float32)
    sample_weights = mask.astype(np.float32)
    donor_cols = np.arange(y.shape[1]) != d6_col
    sample_weights[:, donor_cols] *= arm["donor_weight"]
    active = sample_weights.sum(axis=1) > 0
    if cfg.get("pca_medium"):
        score_dim = cfg["score_dim"]
        score_scaler = StandardScaler().fit(x[train_idx, :score_dim])
        medium_scaler = StandardScaler().fit(x[train_idx, score_dim:])
        pca = PCA(n_components=cfg["pca_medium"], random_state=seed).fit(
            medium_scaler.transform(x[train_idx, score_dim:])
        )
        train_scores = score_scaler.transform(x[train_idx][active, :score_dim])
        held_scores = score_scaler.transform(x[held_idx, :score_dim])
        train_medium = pca.transform(medium_scaler.transform(x[train_idx][active, score_dim:]))
        held_medium = pca.transform(medium_scaler.transform(x[held_idx, score_dim:]))
        pca_scaler = StandardScaler().fit(
            pca.transform(medium_scaler.transform(x[train_idx, score_dim:]))
        )
        train_x = np.column_stack([train_scores, pca_scaler.transform(train_medium)]).astype(
            np.float32
        )
        held_x = np.column_stack([held_scores, pca_scaler.transform(held_medium)]).astype(
            np.float32
        )
    else:
        scaler = StandardScaler().fit(x[train_idx])
        train_x = scaler.transform(x[train_idx][active]).astype(np.float32)
        held_x = scaler.transform(x[held_idx]).astype(np.float32)
    train_targets = normalized[active]
    train_weights = sample_weights[active]
    model = make_model(train_x.shape[1], y.shape[1], cfg).to(device)
    opt = torch.optim.AdamW(
        model.parameters(), lr=cfg["learning_rate"], weight_decay=cfg["weight_decay"]
    )
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=cfg["epochs"])
    rng = np.random.default_rng(seed)
    tensor_x = torch.from_numpy(train_x)
    tensor_y = torch.from_numpy(train_targets)
    tensor_w = torch.from_numpy(train_weights)
    held_tensor = torch.from_numpy(held_x).to(device)
    prediction_epochs = set(cfg["prediction_epochs"])
    held_predictions = []
    final_loss = None
    for epoch in range(1, cfg["epochs"] + 1):
        model.train()
        order = rng.permutation(len(train_x))
        losses = []
        for begin in range(0, len(order), cfg["batch_size"]):
            batch = order[begin : begin + cfg["batch_size"]]
            xb = tensor_x[batch].to(device)
            yb = tensor_y[batch].to(device)
            wb = tensor_w[batch].to(device)
            opt.zero_grad(set_to_none=True)
            pred = model(xb)
            raw = F.smooth_l1_loss(pred, yb, reduction="none", beta=cfg["huber_beta"])
            loss = (raw * wb).sum() / wb.sum()
            loss.backward()
            opt.step()
            losses.append(loss.item())
        scheduler.step()
        final_loss = float(np.mean(losses))
        if epoch in prediction_epochs:
            model.eval()
            with torch.no_grad():
                z = model(held_tensor)[:, d6_col].cpu().numpy()
            held_predictions.append(z * stds[d6_col] + means[d6_col])
    if len(held_predictions) != len(prediction_epochs):
        raise ValueError("Incomplete checkpoints")
    return np.mean(held_predictions, axis=0), {
        "train_compounds": int(active.sum()),
        "d6_train": int(mask[:, d6_col].sum()),
        "donor_train": int(mask[:, donor_cols].sum()),
        "final_train_loss": final_loss,
    }


def main() -> None:
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    torch.set_num_threads(4)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    names = [cfg["target"], *cfg["donors"]]
    d6_col = names.index(cfg["target"])
    direct = (
        pl.read_parquet(inp / "measurements.parquet")
        .filter(pl.col("endpoint").is_in([f"{name}_pIC50_direct_inhibition" for name in names]))
        .with_columns(pl.col("endpoint").str.extract(r"^(CYP[^_]+)", 1).alias("name"))
        .select("compound_id", "name", "value", "conf_low", "conf_high")
    )
    y_table = direct.select("compound_id", "name", "value").pivot(
        on="name", index="compound_id", values="value"
    )
    primary = (
        pl.read_parquet(inp / "regression-labels.parquet")
        .filter(pl.col("head").is_in([f"{name}@0.0000495049505" for name in names]))
        .with_columns(pl.col("head").str.extract(r"^(CYP[^@]+)", 1).alias("name"))
        .select("compound_id", "name", "value")
        .pivot(on="name", index="compound_id", values="value")
        .rename({name: f"{name}_primary" for name in names})
    )
    scores = pl.read_parquet(inp / "scores.parquet")
    medium = pl.read_parquet(inp / "medium.parquet")
    score_cols = sorted(c for c in scores.columns if c.startswith("score_"))
    medium_cols = sorted(c for c in medium.columns if c.startswith("emb_"))
    cols = score_cols + medium_cols
    if len(cols) != 647 or direct.filter(pl.col("name") == cfg["target"]).height != 1493:
        raise ValueError("Unexpected source coverage")
    features = scores.select("compound_id", *score_cols).join(
        medium.select("compound_id", *medium_cols), on="compound_id", validate="1:1"
    )
    output = []
    inventory = []
    for split in cfg["splits"]:
        assignments = pl.read_parquet(inp / f"{split}-assignments.parquet").select(
            "compound_id", "fold", "is_anchor", "group"
        )
        data = (
            features.join(assignments, on="compound_id", validate="1:1")
            .join(y_table, on="compound_id", how="left", validate="1:1")
            .join(primary, on="compound_id", how="left", validate="1:1")
            .sort("compound_id")
        )
        x = data.select(cols).to_numpy().astype(np.float32)
        y = data.select(names).to_numpy().astype(np.float32)
        primary_y = data.select([f"{name}_primary" for name in names]).to_numpy().astype(np.float32)
        if not np.isfinite(x).all() or np.isfinite(y[:, d6_col]).sum() != 1493:
            raise ValueError(f"Invalid {split} feature or label coverage")
        fold_ids = data["fold"].to_numpy()
        anchors = data["is_anchor"].to_numpy()
        for fold in range(5):
            train = np.flatnonzero((fold_ids != fold) | anchors)
            held = np.flatnonzero((fold_ids == fold) & ~anchors & np.isfinite(y[:, d6_col]))
            if len(held) < 100:
                raise ValueError(f"No D6 held rows in {split}/{fold}")
            held_ids = data[held].select("compound_id", "fold", "group")
            truth = direct.filter(pl.col("name") == cfg["target"]).select(
                "compound_id", pl.col("value").alias("y_true"), "conf_low", "conf_high"
            )
            held_frame = held_ids.join(truth, on="compound_id", validate="1:1")
            for arm in cfg["arms"]:
                pred, fit_info = fit_predict(
                    x,
                    y,
                    primary_y,
                    train,
                    held,
                    arm,
                    cfg["model"],
                    cfg["primary_weak_cutoff"],
                    d6_col,
                    cfg["model"]["seed"] + fold,
                    device,
                )
                output.append(
                    held_frame.with_columns(
                        pl.lit(cfg["target"] + "_pIC50_direct_inhibition").alias("endpoint"),
                        pl.lit(split).alias("split_id"),
                        pl.lit(arm["name"]).alias("arm"),
                        pl.Series("prediction", pred),
                    )
                )
                inventory.append({"split_id": split, "fold": fold, "arm": arm["name"], **fit_info})
                print(split, fold, arm["name"], fit_info, flush=True)
    oof = pl.concat(output).sort("split_id", "arm", "compound_id")
    expected = len(cfg["splits"]) * len(cfg["arms"]) * 1457
    if oof.height != expected or oof.select("split_id", "arm", "compound_id").is_duplicated().any():
        raise ValueError("Incomplete D6 OOF")
    rows = []
    for (split, arm), frame in oof.group_by("split_id", "arm"):
        y_true = frame["y_true"].to_numpy()
        pred = frame["prediction"].to_numpy()
        rows.append(
            {
                "split_id": split,
                "arm": arm,
                "n": frame.height,
                "ST_RAE": st_rae(frame),
                "MAE": float(np.abs(y_true - pred).mean()),
                "Spearman": float(spearmanr(y_true, pred).statistic),
                "mean_prediction": float(pred.mean()),
            }
        )
    oof.write_parquet(out / "oof.parquet")
    pl.DataFrame(rows).sort("split_id", "arm").write_csv(out / "metrics.csv")
    pl.DataFrame(inventory).sort("split_id", "fold", "arm").write_csv(out / "inventory.csv")
    (out / "summary.json").write_text(
        json.dumps({"config": cfg, "device": str(device)}, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
