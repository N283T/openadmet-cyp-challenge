"""Fresh D6 pilot with paired-measurement learning and nested operating points.

No historical score/model/call-count input is accepted. Auxiliary paired measurements
are restricted to training molecules at every outer and inner fitting boundary.
"""

import json
import os
from pathlib import Path

import joblib
import numpy as np
import polars as pl
from lightgbm import LGBMClassifier, LGBMRegressor
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold
from sklearn.naive_bayes import ComplementNB
from sklearn.svm import SVC


def bootstrap_mcc(y, score, thresholds, samples, seed):
    """Evaluate all thresholds with the same seed-zero empirical bootstrap draws."""
    order = np.argsort(score, kind="stable")
    y = np.asarray(y, dtype=int)[order]
    cuts = np.searchsorted(np.asarray(score)[order], thresholds, side="right")
    draws = np.random.default_rng(seed).integers(0, len(y), (samples, len(y)))
    weights = np.zeros((samples, len(y)), dtype=np.int32)
    np.add.at(weights, (np.arange(samples)[:, None], draws), 1)
    cp = np.pad(np.cumsum(weights * y, axis=1), ((0, 0), (1, 0)))
    cn = np.pad(np.cumsum(weights * (1 - y), axis=1), ((0, 0), (1, 0)))
    fn, tn = cp[:, cuts], cn[:, cuts]
    tp, fp = cp[:, -1, None] - fn, cn[:, -1, None] - tn
    denom = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    values = np.divide(tp * tn - fp * fn, denom, out=np.zeros_like(denom), where=denom > 0)
    return values.mean(axis=0), values.std(axis=0, ddof=1)


def folds(indices, y, groups, strategy, count, seed):
    # Missing endpoint labels are a distinct stratum, never labelled negatives.
    strata = np.where(np.isfinite(y[indices]), y[indices], 2).astype(int)
    if strategy in {"umap", "scaffold", "similarity"}:
        splitter = StratifiedGroupKFold(n_splits=count, shuffle=True, random_state=seed)
        iterator = splitter.split(indices, strata, groups[indices])
    else:
        splitter = StratifiedKFold(n_splits=count, shuffle=True, random_state=seed)
        iterator = splitter.split(indices, strata)
    return [(indices[a], indices[b]) for a, b in iterator]


def fit_multiassay(train, xfp, tasks, seed, cfg):
    """Fresh masked multitask molecular network; every scale is fitted inside train."""
    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
    import torch
    from torch import nn

    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.use_deterministic_algorithms(True)
    torch.set_float32_matmul_precision("high")
    device = "cuda" if torch.cuda.is_available() else "cpu"
    target = tasks[train].copy()
    mask = np.isfinite(target)
    class_weights = np.ones(len(train), dtype=np.float32)
    if cfg.get("selection_weighting"):
        from cyp_primary_supervision import selection_weights

        class_weights = selection_weights(train, tasks, cfg.get("propensity_weight_cap", 5.0))
    means = np.nanmean(target[:, 4:], axis=0)
    scales = np.maximum(np.nanstd(target[:, 4:], axis=0), 0.1)
    target[:, 4:] = (target[:, 4:] - means) / scales
    target[~mask] = 0
    inputs = torch.as_tensor(xfp[train], device=device)
    target = torch.as_tensor(target, dtype=torch.float32, device=device)
    mask = torch.as_tensor(mask, dtype=torch.float32, device=device)
    denominator = mask.sum(dim=0).clamp_min(1)
    class_weights = torch.as_tensor(class_weights, device=device)
    network = nn.Sequential(
        nn.Linear(2048, 256),
        nn.ReLU(),
        nn.Dropout(0.15),
        nn.Linear(256, 64),
        nn.ReLU(),
        nn.Linear(64, tasks.shape[1]),
    ).to(device)
    optimizer = torch.optim.AdamW(network.parameters(), lr=0.002, weight_decay=0.001)
    for _ in range(cfg.get("network_epochs", 100)):
        network.train()
        optimizer.zero_grad(set_to_none=True)
        prediction = network(inputs)
        binary = nn.functional.binary_cross_entropy_with_logits(
            prediction[:, :4], target[:, :4], reduction="none"
        )
        regression = (prediction[:, 4:] - target[:, 4:]).square()
        loss_binary = (binary * mask[:, :4]).sum(dim=0) / denominator[:4]
        loss_aux = (regression * mask[:, 4:]).sum(dim=0) / denominator[4:]
        auxiliary = loss_aux[cfg.get("auxiliary_indices", list(range(len(loss_aux))))].mean()
        loss = (
            (binary[:, 0] * mask[:, 0] * class_weights).sum()
            / (mask[:, 0] * class_weights).sum().clamp_min(1)
            + cfg.get("other_enzyme_weight", 0.3) * loss_binary[1:].mean()
            + cfg.get("continuous_aux_weight", 0.2) * auxiliary
        )
        loss.backward()
        optimizer.step()
    network.eval()
    return dict(
        network=network.cpu(),
        means=means,
        scales=scales,
        seed=seed,
        device_used=device,
        epochs=cfg.get("network_epochs", 100),
    )


def fit_model(kind, train, xfp, x, y, tdi, shift, seed, weights, subtype, tasks, cfg):
    if kind.startswith("chemeleon_"):
        embedding = x[:, -cfg["embedding_dimensions"] :]
        mean = embedding[train].mean(axis=0)
        scale = np.maximum(embedding[train].std(axis=0), 0.01)
        normalized = ((embedding - mean) / scale).astype(np.float32)
        network_tasks = tasks.copy()
        network_tasks[:, 0] = y
        network_cfg = {
            **cfg,
            "other_enzyme_weight": 0.0,
            "continuous_aux_weight": 0.2 if "primary" in kind else 0.0,
            "auxiliary_indices": [12, 13, 14, 15],
        }
        model = fit_multiassay(train, normalized, network_tasks, seed, network_cfg)
        model.update(feature_mean=mean, feature_scale=scale)
        return model
    if kind in {"potency_all", "potency_conditional", "potency_gate"}:
        from cyp_primary_supervision import fit_potency_gate

        return fit_potency_gate(kind, train, xfp, y, tasks, seed, cfg)
    if kind in {"activity_all", "activity_selected", "gate_selected"}:
        from cyp_primary_supervision import fit_activity_gate

        return fit_activity_gate(train, xfp, y, tasks, seed, cfg)
    if kind in {"primary_d6_network", "primary_cross_cyp_network", "primary_selection_network"}:
        network_cfg = {
            **cfg,
            "other_enzyme_weight": 0.0,
            "continuous_aux_weight": 0.2,
            "auxiliary_indices": (
                [12] if kind != "primary_cross_cyp_network" else [12, 13, 14, 15]
            ),
            "selection_weighting": kind == "primary_selection_network",
        }
        return fit_multiassay(train, xfp, tasks, seed, network_cfg)
    if kind in {"fresh_multiassay_network", "fresh_d6_network"}:
        network_cfg = (
            cfg
            if kind == "fresh_multiassay_network"
            else {**cfg, "other_enzyme_weight": 0.0, "continuous_aux_weight": 0.0}
        )
        return fit_multiassay(train, xfp, tasks, seed, network_cfg)
    if kind in {
        "historical_morgan_lr_c005",
        "historical_morgan_rbf_c1",
        "historical_morgan_complement_a1",
    }:
        rows = train[np.isfinite(y[train])]
        factory = {
            "historical_morgan_lr_c005": lambda: LogisticRegression(
                C=0.05, max_iter=4000, solver="liblinear", random_state=seed
            ),
            "historical_morgan_rbf_c1": lambda: SVC(C=1.0, kernel="rbf"),
            "historical_morgan_complement_a1": lambda: ComplementNB(alpha=1.0),
        }
        model = factory[kind]()
        model.fit(xfp[rows], y[rows])
        return model
    params = dict(
        n_estimators=240,
        num_leaves=15,
        min_child_samples=35,
        learning_rate=0.035,
        reg_lambda=5.0,
        verbosity=-1,
        n_jobs=4,
        random_state=seed,
        deterministic=True,
        force_col_wise=True,
    )
    if kind in {"morgan_logistic", "confidence_weighted_logistic"}:
        rows = train[np.isfinite(y[train])]
        model = LogisticRegression(C=1.0, max_iter=1200, solver="liblinear", random_state=seed)
        model.fit(
            xfp[rows],
            y[rows],
            sample_weight=(
                weights[rows] / weights[rows].mean()
                if kind == "confidence_weighted_logistic"
                else None
            ),
        )
    elif kind in {"regime_multiclass", "confidence_weighted_regime"}:
        rows = train[np.isfinite(y[train])]
        model = LogisticRegression(C=1.0, max_iter=1200, solver="lbfgs", random_state=seed)
        model.fit(
            xfp[rows],
            subtype[rows],
            sample_weight=(
                weights[rows] / weights[rows].mean()
                if kind == "confidence_weighted_regime"
                else None
            ),
        )
    elif kind == "morgan_rdkit_classifier":
        rows = train[np.isfinite(y[train])]
        model = LGBMClassifier(**params)
        model.fit(x[rows], y[rows])
    elif kind == "paired_measurement_margin":
        model = []
        for target in [tdi, shift]:
            rows = train[np.isfinite(target[train])]
            regressor = LGBMRegressor(**params)
            regressor.fit(x[rows], target[rows])
            model.append(regressor)
    else:
        raise ValueError(f"Unknown model: {kind}")
    return model


def predict(kind, model, rows, xfp, x):
    if kind.startswith("chemeleon_"):
        import torch

        dimensions = len(model["feature_mean"])
        features = (
            (x[rows, -dimensions:] - model["feature_mean"]) / model["feature_scale"]
        ).astype(np.float32)
        with torch.no_grad():
            return torch.sigmoid(model["network"](torch.as_tensor(features))[:, 0]).numpy()
    if kind in {
        "activity_all",
        "activity_selected",
        "gate_selected",
        "potency_all",
        "potency_conditional",
        "potency_gate",
    }:
        import torch

        with torch.no_grad():
            probabilities = torch.sigmoid(model["network"](torch.as_tensor(xfp[rows])))
            return (
                probabilities[:, 0] * probabilities[:, 1]
                if kind in {"gate_selected", "potency_gate"}
                else probabilities[:, 0]
            ).numpy()
    if kind == "historical_morgan_rbf_c1":
        return model.decision_function(xfp[rows])
    if kind in {"historical_morgan_lr_c005", "historical_morgan_complement_a1"}:
        return model.predict_proba(xfp[rows])[:, 1]
    if kind in {
        "fresh_multiassay_network",
        "fresh_d6_network",
        "primary_d6_network",
        "primary_cross_cyp_network",
        "primary_selection_network",
    }:
        import torch

        with torch.no_grad():
            return torch.sigmoid(model["network"](torch.as_tensor(xfp[rows]))[:, 0]).numpy()
    if kind in {"morgan_logistic", "confidence_weighted_logistic"}:
        return model.predict_proba(xfp[rows])[:, 1]
    if kind in {"regime_multiclass", "confidence_weighted_regime"}:
        return model.predict_proba(xfp[rows])[:, model.classes_ > 0].sum(axis=1)
    if kind == "morgan_rdkit_classifier":
        return model.predict_proba(x[rows])[:, 1]
    potency = model[0].predict(x[rows])
    shift = model[1].predict(x[rows])
    return np.minimum(potency - 4 - np.log10(2), shift - np.log10(2))


def metrics(y, score, calls, samples):
    keep = np.isfinite(y)
    y, score, calls = y[keep].astype(int), score[keep], calls[keep]
    mean, std = bootstrap_mcc(y, calls.astype(float), [0.5], samples, 0)
    return dict(
        n=len(y),
        positives=int(y.sum()),
        predicted_positive=int(calls.sum()),
        ap=float(average_precision_score(y, score)),
        auc=float(roc_auc_score(y, score)) if len(np.unique(y)) > 1 else None,
        mcc=float(matthews_corrcoef(y, calls)),
        bootstrap_mcc=float(mean[0]),
        bootstrap_std=float(std[0]),
    )


def main():
    RDLogger.DisableLog("rdApp.*")
    source, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    endpoint, seed = cfg["endpoint"], cfg["seed"]
    target_prefix = endpoint.removesuffix("_is_TDI")
    compounds = pl.read_parquet(source / "compounds.parquet").sort("compound_id")
    aliases = pl.read_parquet(source / "aliases.parquet")
    names = aliases.group_by("compound_id").agg(pl.col("molecule_name").sort().first())
    test_ids = aliases.filter(pl.col("role") == "test")["compound_id"].unique().to_list()
    measured = pl.read_parquet(source / "measurements.parquet")
    assert measured.select("compound_id", "endpoint").is_duplicated().sum() == 0
    wide = measured.pivot(on="endpoint", index="compound_id", values="value")
    table = compounds.join(names, on="compound_id", how="left").join(
        wide, on="compound_id", how="left"
    )
    y = table[endpoint].to_numpy()
    direct = table[f"{target_prefix}_pIC50_direct_inhibition"].to_numpy()
    tdi = table[f"{target_prefix}_pIC50_TDI_condition"].to_numpy()
    shift = tdi - direct
    blind = np.flatnonzero(table["compound_id"].is_in(test_ids).to_numpy())
    development = np.flatnonzero(~table["compound_id"].is_in(test_ids).to_numpy())
    assert len(blind) == 750 and not np.isfinite(y[blind]).any()
    regimes = np.full(len(table), "no_D6_label", dtype=object)
    paired = np.isfinite(direct) & np.isfinite(tdi)
    labelled = np.isfinite(y)
    regimes[labelled] = "unresolved_or_unpaired"
    regimes[labelled & paired & (direct > 4)] = "measured_direct_active"
    regimes[labelled & paired & (direct < 4) & (tdi < 4)] = "assigned_inactive"
    regimes[labelled & paired & (direct < 4) & (tdi > 4 + np.log10(2))] = "inferred_positive"
    bounds = []
    for condition, name in [("direct_inhibition", "direct"), ("TDI_condition", "tdi")]:
        interval = measured.filter(
            pl.col("endpoint") == f"{target_prefix}_pIC50_{condition}"
        ).select(
            "compound_id",
            pl.col("conf_low").alias(f"{name}_low"),
            pl.col("conf_high").alias(f"{name}_high"),
        )
        ordered = table.select("compound_id").join(interval, on="compound_id", how="left")
        bounds.extend([ordered[f"{name}_low"].to_numpy(), ordered[f"{name}_high"].to_numpy()])
    direct_low, direct_high, tdi_low, tdi_high = bounds
    strong_shift_positive = tdi_low - direct_high > np.log10(2)
    strong_shift_negative = tdi_high - direct_low <= np.log10(2)
    reliable = (
        (
            (regimes == "measured_direct_active")
            & (((y == 1) & strong_shift_positive) | ((y == 0) & strong_shift_negative))
        )
        | ((regimes == "inferred_positive") & (tdi_low > 4 + np.log10(2)))
        | ((regimes == "assigned_inactive") & (direct_high < 4) & (tdi_high < 4))
    )
    weights = np.where(reliable, 1.0, cfg.get("ambiguous_weight", 0.25))
    subtype = np.where(y == 1, np.where(regimes == "inferred_positive", 2, 1), 0)
    enzymes = [target_prefix] + [
        e for e in ["CYP2D6", "CYP3A4", "CYP1A2", "CYP2C9"] if e != target_prefix
    ]
    tasks = [table[f"{enzyme}_is_TDI"].to_numpy() for enzyme in enzymes]
    for enzyme in enzymes:
        active_direct = np.maximum(table[f"{enzyme}_pIC50_direct_inhibition"].to_numpy(), 4)
        active_tdi = np.maximum(table[f"{enzyme}_pIC50_TDI_condition"].to_numpy(), 4)
        tasks.extend([active_direct, active_tdi, active_tdi - active_direct])
    tasks = np.stack(tasks, axis=1)
    if cfg.get("primary_screen"):
        from cyp_primary_supervision import load_primary

        primary, primary_rows = load_primary(source / cfg["primary_screen"], aliases, table)
        assert not np.isfinite(primary[blind]).any()
        primary_rows.write_parquet(out / "primary-source-audit.parquet")
        tasks = np.concatenate([tasks, primary], axis=1)
    audit = table.select("compound_id", "molecule_name").with_columns(
        pl.Series("label", y),
        pl.Series("direct_pIC50", direct),
        pl.Series("tdi_pIC50", tdi),
        pl.Series("shift", shift),
        pl.Series("regime", regimes),
        pl.Series("training_weight", weights),
        pl.Series("positive_subtype", subtype),
        pl.Series("is_blind", np.isin(np.arange(len(table)), blind)),
    )
    audit.write_parquet(out / "source-audit.parquet")
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    # 3D descriptors are not requested; the registered list here is RDKit's 2D set.
    descriptor_names = [name for name, _ in Descriptors._descList]
    fp = np.zeros((len(table), 2048), dtype=np.float32)
    desc = np.full((len(table), len(descriptor_names)), np.nan, dtype=np.float32)
    groups = []
    fingerprints = []
    for i, smiles in enumerate(table["standardized_smiles"]):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Unparseable canonical molecule: {table['compound_id'][i]}")
        fingerprint = generator.GetFingerprint(mol)
        fingerprints.append(fingerprint)
        DataStructs.ConvertToNumpyArray(fingerprint, fp[i])
        for j, (_, function) in enumerate(Descriptors._descList):
            try:
                desc[i, j] = function(mol)
            except (ValueError, RuntimeError, OverflowError, ZeroDivisionError):
                pass
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
        # Acyclic molecules have no common scaffold; avoid one giant empty group.
        groups.append(scaffold or f"acyclic:{table['compound_id'][i]}")
    groups = np.asarray(groups)
    desc[~np.isfinite(desc)] = np.nan
    x = np.concatenate([fp, desc], axis=1)
    if cfg.get("embedding_file"):
        embeddings = pl.read_parquet(source / cfg["embedding_file"])
        columns = [c for c in embeddings.columns if c.startswith("emb_")]
        assert len(columns) == cfg["embedding_dimensions"] == 2048
        aligned = (
            table.select("compound_id")
            .join(embeddings, on="compound_id", how="left", validate="1:1", maintain_order="left")
            .select(columns)
            .to_numpy()
            .astype(np.float32)
        )
        assert np.isfinite(aligned).all() and embeddings.height == len(table)
        x = np.concatenate([x, aligned], axis=1)
    if "similarity" in cfg["split_strategies"]:
        from scipy.sparse import csr_matrix
        from scipy.sparse.csgraph import connected_components

        edge_rows, edge_columns = [], []
        train_fingerprints = [fingerprints[i] for i in development]
        for i, fingerprint in enumerate(train_fingerprints):
            similarity = np.asarray(
                DataStructs.BulkTanimotoSimilarity(fingerprint, train_fingerprints[:i])
            )
            neighbors = np.flatnonzero(similarity >= cfg["similarity_threshold"])
            edge_rows.extend([i] * len(neighbors))
            edge_columns.extend(neighbors.tolist())
        graph = csr_matrix(
            (np.ones(len(edge_rows)), (edge_rows, edge_columns)),
            shape=(len(development), len(development)),
        )
        count, cluster = connected_components(graph, directed=False)
        groups = np.full(len(table), "blind", dtype=object)
        groups[development] = [f"similarity:{i}" for i in cluster]
        cluster_table = (
            table[development]
            .select("compound_id", "molecule_name")
            .with_columns(
                pl.Series("group", groups[development]), pl.Series("label", y[development])
            )
        )
        cluster_table.write_parquet(out / "chemical-groups.parquet")
        counts = cluster_table.group_by("group").agg(
            pl.len().alias("molecules"), pl.col("label").is_finite().sum().alias("D6_labels")
        )
        counts.write_csv(out / "chemical-group-sizes.csv")
        info = dict(
            threshold=cfg["similarity_threshold"],
            groups=count,
            edges=len(edge_rows),
            largest_group=int(counts["molecules"].max()),
            labelled_groups=int((counts["D6_labels"] > 0).sum()),
            guarantee="No train-validation pair has Morgan Tanimoto >= threshold",
        )
        (out / "chemical-group-summary.json").write_text(json.dumps(info, indent=2) + "\n")
        print(json.dumps(info), flush=True)
    (out / "representation.json").write_text(
        json.dumps(
            dict(
                descriptors=descriptor_names,
                fingerprint_bits=2048,
                radius=2,
                no_pretrained_assets=not bool(cfg.get("embedding_file")),
                pretrained_features=cfg.get("embedding_file"),
            ),
            indent=2,
        )
        + "\n"
    )
    if cfg.get("deploy_nested_policy"):
        from deploy_cyp_chemeleon_candidate import deploy

        deploy(table, fp, x, y, tdi, shift, weights, subtype, tasks, development, blind, cfg)
        return
    if cfg.get("nested_bag_validation"):
        from validate_cyp_d6_bag_policy import evaluate

        evaluate(table, fp, x, y, tdi, shift, weights, subtype, tasks, development, blind, cfg)
        return
    results, predictions, operating_points, assignments = [], [], [], []
    blind_fold_predictions = []
    for strategy in cfg["split_strategies"]:
        if cfg.get("fixed_assignments"):
            saved = pl.read_parquet(source / cfg["fixed_assignments"]).filter(
                pl.col("strategy") == strategy
            )
            ordered = table.select("compound_id").join(
                saved.select("compound_id", "fold", "group"),
                on="compound_id",
                how="left",
                maintain_order="left",
            )
            membership = ordered["fold"].to_numpy()
            groups = ordered["group"].fill_null("blind").to_numpy()
            assert np.isfinite(membership[development]).all()
            outer = [
                (
                    development[membership[development] != f],
                    development[membership[development] == f],
                )
                for f in range(cfg["outer_folds"])
            ]
        else:
            outer = folds(development, y, groups, strategy, cfg["outer_folds"], seed)
        for fold, (train, test) in enumerate(outer):
            if strategy in {"umap", "scaffold", "similarity"}:
                assert set(groups[train]).isdisjoint(groups[test])
            for row in test:
                assignments.append(
                    dict(
                        compound_id=table["compound_id"][int(row)],
                        split_id=strategy,
                        fold=fold,
                        group=groups[row],
                    )
                )
            inner = folds(train, y, groups, strategy, cfg["inner_folds"], seed + 100 + fold)
            for kind in cfg["models"]:
                inner_score = np.full(len(table), np.nan)
                for number, (fit, valid) in enumerate(inner):
                    model = fit_model(
                        kind,
                        fit,
                        fp,
                        x,
                        y,
                        tdi,
                        shift,
                        seed + fold * 10 + number,
                        weights,
                        subtype,
                        tasks,
                        cfg,
                    )
                    inner_score[valid] = predict(kind, model, valid, fp, x)
                observed = train[np.isfinite(y[train])]
                grid = (
                    np.linspace(-1.0, 1.0, 101)
                    if kind == "paired_measurement_margin"
                    else np.linspace(0.01, 0.99, 99)
                )
                grid = np.r_[-np.inf, grid, np.inf]
                if cfg.get("threshold_grid") == "inner_quantiles":
                    grid = np.r_[
                        -np.inf,
                        np.unique(np.quantile(inner_score[observed], np.linspace(0, 1, 301))),
                        np.inf,
                    ]
                mean, std = bootstrap_mcc(
                    y[observed], inner_score[observed], grid, cfg["bootstrap_samples"], 0
                )
                best = int(np.argmax(mean))
                threshold = float(grid[best])
                model = fit_model(
                    kind, train, fp, x, y, tdi, shift, seed + fold, weights, subtype, tasks, cfg
                )
                score = predict(kind, model, test, fp, x)
                if cfg.get("save_blind_folds"):
                    blind_score = predict(kind, model, blind, fp, x)
                    for row, s in zip(blind, blind_score, strict=True):
                        blind_fold_predictions.append(
                            dict(
                                compound_id=table["compound_id"][int(row)],
                                endpoint=endpoint,
                                split_id=strategy,
                                fold=fold,
                                model=kind,
                                score=float(s),
                            )
                        )
                calls = score > threshold
                result = metrics(y[test], score, calls, cfg["bootstrap_samples"])
                results.append(dict(split_id=strategy, fold=fold, model=kind, **result))
                operating_points.append(
                    dict(
                        split_id=strategy,
                        fold=fold,
                        model=kind,
                        threshold=threshold,
                        inner_bootstrap_mcc=float(mean[best]),
                        inner_bootstrap_std=float(std[best]),
                        fitting_molecules=len(train),
                    )
                )
                for row, s, call in zip(test, score, calls, strict=True):
                    predictions.append(
                        dict(
                            compound_id=table["compound_id"][int(row)],
                            Molecule_Name=table["molecule_name"][int(row)],
                            endpoint=endpoint,
                            split_id=strategy,
                            fold=fold,
                            model=kind,
                            score=float(s),
                            call=bool(call),
                            label=float(y[row]),
                            regime=regimes[row],
                        )
                    )
                print(json.dumps(results[-1]), flush=True)
        # Each family is fit from zero on all development data. Its threshold is
        # selected on fresh outer OOF scores, so blind output has no fitted call budget.
        for kind in cfg["models"]:
            frame = pl.DataFrame(predictions).filter(
                (pl.col("split_id") == strategy) & (pl.col("model") == kind)
            )
            observed = frame.filter(pl.col("label").is_finite())
            score, labels, calls = (
                observed["score"].to_numpy(),
                observed["label"].to_numpy(),
                observed["call"].to_numpy(),
            )
            results.append(
                dict(
                    split_id=strategy,
                    fold=-1,
                    model=kind,
                    **metrics(labels, score, calls, cfg["bootstrap_samples"]),
                )
            )
            grid = (
                np.linspace(-1, 1, 101)
                if kind == "paired_measurement_margin"
                else np.linspace(0.01, 0.99, 99)
            )
            grid = np.r_[-np.inf, grid, np.inf]
            if cfg.get("threshold_grid") == "inner_quantiles":
                grid = np.r_[-np.inf, np.unique(np.quantile(score, np.linspace(0, 1, 301))), np.inf]
            means, _ = bootstrap_mcc(labels, score, grid, cfg["bootstrap_samples"], 0)
            threshold = float(grid[int(np.argmax(means))])
            model = fit_model(
                kind, development, fp, x, y, tdi, shift, seed, weights, subtype, tasks, cfg
            )
            scores = predict(kind, model, blind, fp, x)
            blind_table = (
                table[blind]
                .select("compound_id", "molecule_name")
                .with_columns(
                    pl.lit(endpoint).alias("endpoint"),
                    pl.lit("blind").alias("split_id"),
                    pl.lit(-1).alias("fold"),
                    pl.Series("score", scores),
                    pl.Series("call", scores > threshold),
                )
            )
            blind_table.write_parquet(out / f"blind-{strategy}-{kind}.parquet")
            joblib.dump(model, out / f"model-{strategy}-{kind}.joblib")
            operating_points.append(
                dict(
                    split_id=strategy,
                    fold=-1,
                    model=kind,
                    threshold=threshold,
                    blind_positive=int((scores > threshold).sum()),
                    note="Full-development threshold; pooled OOF fit is an optimistic MCC estimate",
                )
            )
    pl.DataFrame(predictions).write_parquet(out / "oof.parquet")
    if blind_fold_predictions:
        pl.DataFrame(blind_fold_predictions).write_parquet(out / "blind-folds.parquet")
    pl.DataFrame(assignments).write_parquet(out / "assignments.parquet")
    pl.DataFrame(results).write_csv(out / "metrics.csv")
    pl.DataFrame(operating_points, infer_schema_length=None).write_csv(out / "operating-points.csv")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                config=cfg,
                aggregate=[row for row in results if row["fold"] == -1],
                operating_points=[row for row in operating_points if row["fold"] == -1],
                caveats=[
                    "Scored blind label mask unknown",
                    "Pilot cannot establish transfer to Live",
                    "Regression point estimates below pIC50 4 are censored/noisy",
                    "Scaffold split leaves acyclic compounds as singletons",
                    "Bootstrap resamples use NumPy generator rather than official sklearn RNG",
                ],
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
