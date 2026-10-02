"""Primary-screen supervision and training-only label-selection propensity."""

import numpy as np
import polars as pl
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures, StandardScaler


def load_primary(path, aliases, table):
    mapping = aliases.select("molecule_name", "compound_id").unique()
    assert not mapping["molecule_name"].is_duplicated().any()
    raw = pl.read_csv(path).rename({"Molecule_Name": "molecule_name"})
    linked = raw.join(mapping, on="molecule_name", how="left", validate="m:1")
    if linked["compound_id"].null_count():
        raise ValueError("Primary-screen compounds missing canonical identity")
    concentrations = linked["concentration_M"].unique().to_list()
    assert len(concentrations) == 1 and 0.000049 <= concentrations[0] <= 0.000051
    grouped = linked.group_by("compound_id", "enzyme").agg(
        pl.col("log2fc_estimate").median().alias("log2fc"),
        pl.len().alias("source_rows"),
    )
    wide = grouped.pivot(on="enzyme", index="compound_id", values="log2fc")
    ordered = table.select("compound_id").join(
        wide, on="compound_id", how="left", maintain_order="left"
    )
    matrix = np.column_stack(
        [ordered[e].to_numpy() for e in ["CYP2D6", "CYP3A4", "CYP1A2", "CYP2C9"]]
    )
    return matrix, grouped


def propensity(fit_primary, fit_available, predict_primary):
    observed = np.isfinite(fit_primary)
    median = float(np.median(fit_primary[observed])) if observed.any() else 0.0

    def features(values):
        return np.column_stack(
            [np.clip(np.where(np.isfinite(values), values, median), -8, 2), ~np.isfinite(values)]
        ).astype(float)

    model = make_pipeline(
        StandardScaler(),
        PolynomialFeatures(3, include_bias=False),
        LogisticRegression(C=1.0, max_iter=1000, solver="lbfgs"),
    )
    model.fit(features(fit_primary), np.asarray(fit_available, dtype=int))
    return model.predict_proba(features(predict_primary))[:, 1]


def selection_weights(train, tasks, cap=5.0):
    """Weights correct selection only under a missing-at-random assumption, not known truth."""
    available = np.isfinite(tasks[train, 0])
    primary = tasks[train, 16]
    probability = propensity(primary, available, primary)
    weights = 1 / np.maximum(probability, 1 / cap)
    weights /= weights[available].mean()
    return weights.astype(np.float32)


def fit_activity_gate(train, fp, y, tasks, seed, cfg):
    """Learn measured screen activity and observed TDI separately; no missing class imputation."""
    import torch
    from torch import nn

    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.use_deterministic_algorithms(True)
    torch.set_float32_matmul_precision(cfg.get("float32_matmul_precision", "highest"))
    device = "cuda" if torch.cuda.is_available() else "cpu"
    primary = tasks[train, 16 + cfg["primary_enzyme_index"]]
    mask = np.column_stack([np.isfinite(y[train]), np.isfinite(primary)])
    target = np.column_stack(
        [
            np.nan_to_num(y[train]),
            (primary < cfg["primary_activity_cut"]).astype(float),
        ]
    )
    network = nn.Sequential(
        nn.Linear(2048, 256),
        nn.ReLU(),
        nn.Dropout(0.15),
        nn.Linear(256, 64),
        nn.ReLU(),
        nn.Linear(64, 2),
    ).to(device)
    inputs = torch.as_tensor(fp[train], device=device)
    target = torch.as_tensor(target, dtype=torch.float32, device=device)
    mask = torch.as_tensor(mask, dtype=torch.float32, device=device)
    denominator = mask.sum(dim=0).clamp_min(1)
    optimizer = torch.optim.AdamW(network.parameters(), lr=0.002, weight_decay=0.001)
    for _ in range(cfg["network_epochs"]):
        network.train()
        optimizer.zero_grad(set_to_none=True)
        prediction = network(inputs)
        loss = nn.functional.binary_cross_entropy_with_logits(prediction, target, reduction="none")
        losses = (loss * mask).sum(dim=0) / denominator
        (losses[0] + cfg["primary_activity_loss_weight"] * losses[1]).backward()
        optimizer.step()
    network.eval()
    return {"network": network.cpu()}


def fit_potency_gate(kind, train, fp, y, tasks, seed, cfg):
    """Separate a measured TDI-potency necessary condition from conditional TDI class."""
    import torch
    from torch import nn

    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.use_deterministic_algorithms(True)
    torch.set_float32_matmul_precision(cfg.get("float32_matmul_precision", "highest"))
    device = "cuda" if torch.cuda.is_available() else "cpu"
    potency = tasks[train, 5]
    active = potency > 4 + np.log10(2)
    assert not np.any((y[train] == 1) & ~active)
    primary = tasks[train, 16 + cfg["primary_enzyme_index"]]
    primary_mean, primary_scale = np.nanmean(primary), max(np.nanstd(primary), 0.1)
    class_mask = np.isfinite(y[train])
    if kind != "potency_all":
        class_mask &= active
    mask = np.column_stack([class_mask, np.isfinite(potency), np.isfinite(primary)])
    target = np.column_stack(
        [
            np.nan_to_num(y[train]),
            active.astype(float),
            np.nan_to_num((primary - primary_mean) / primary_scale),
        ]
    )
    network = nn.Sequential(
        nn.Linear(2048, 256),
        nn.ReLU(),
        nn.Dropout(0.15),
        nn.Linear(256, 64),
        nn.ReLU(),
        nn.Linear(64, 3),
    ).to(device)
    inputs = torch.as_tensor(fp[train], device=device)
    target = torch.as_tensor(target, dtype=torch.float32, device=device)
    mask = torch.as_tensor(mask, dtype=torch.float32, device=device)
    denominator = mask.sum(dim=0).clamp_min(1)
    optimizer = torch.optim.AdamW(network.parameters(), lr=0.002, weight_decay=0.001)
    for _ in range(cfg["network_epochs"]):
        network.train()
        optimizer.zero_grad(set_to_none=True)
        prediction = network(inputs)
        binary = nn.functional.binary_cross_entropy_with_logits(
            prediction[:, :2], target[:, :2], reduction="none"
        )
        losses = (binary * mask[:, :2]).sum(dim=0) / denominator[:2]
        auxiliary = ((prediction[:, 2] - target[:, 2]).square() * mask[:, 2]).sum() / denominator[2]
        (losses.sum() + 0.2 * auxiliary).backward()
        optimizer.step()
    network.eval()
    return {"network": network.cpu(), "primary_mean": primary_mean, "primary_scale": primary_scale}
