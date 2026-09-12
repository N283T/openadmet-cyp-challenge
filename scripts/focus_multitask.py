"""Focus-weighted fold-safe external600 pIC50 and auxiliary multitask training.

Outer nonanchor compounds lose ALL head labels during fitting and checkpoint
selection. Anchors stay in training. Inner validation separates supplied molecular
groups; groups containing anchors cannot enter inner validation.
"""

import copy
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
from chemprop import featurizers
from chemprop import nn as cnn
from chemprop.data import BatchMolGraph
from rdkit import Chem


def input_path(root, name):
    paths = [path for path in root.iterdir() if path.name.partition("-")[2] == name]
    if len(paths) != 1:
        raise ValueError(f"Expected one input {name}; got {paths}")
    return paths[0]


def encoder_architecture(checkpoint):
    source = checkpoint.get("config", {})
    medium = source.get("architecture") == "medium"
    architecture = dict(
        width=source.get("width", 600),
        depth=source.get("depth", 4 if medium else 5),
        aggregation=source.get("aggregation", "norm" if medium else "mean"),
        aggregation_norm=source.get("aggregation_norm", 100.0) if medium else None,
    )
    assert architecture["width"] == 600
    assert architecture["depth"] == (4 if medium else 5)
    assert architecture["aggregation"] == ("norm" if medium else "mean")
    if medium:
        assert architecture["aggregation_norm"] == 100.0
    return architecture


class MultitaskModel(torch.nn.Module):
    def __init__(self, encoder_state, dropout, n_heads=8, architecture=None):
        super().__init__()
        architecture = architecture or encoder_architecture({})
        self.mp = cnn.BondMessagePassing(d_h=600, depth=architecture["depth"])
        self.mp.load_state_dict(encoder_state, strict=True)
        self.aggregation = (
            cnn.NormAggregation(norm=architecture["aggregation_norm"])
            if architecture["aggregation"] == "norm"
            else cnn.MeanAggregation()
        )
        self.head = torch.nn.Sequential(
            torch.nn.Linear(600, 256),
            torch.nn.ReLU(),
            torch.nn.Dropout(dropout),
            torch.nn.Linear(256, n_heads),
        )

    def embed(self, graph):
        return self.aggregation(self.mp(graph), graph.batch)

    def forward(self, graph):
        return self.head(self.embed(graph))


def masked_loss(prediction, truth, weights, kind="mse"):
    mask = torch.isfinite(truth)
    error = prediction - torch.nan_to_num(truth)
    if kind == "mse":
        error = error.square()
    elif kind == "mae":
        error = error.abs()
    else:
        raise ValueError(f"Unknown training loss: {kind}")
    per_head = (error * mask).sum(0)
    per_head = per_head / mask.sum(0).clamp_min(1)
    present = mask.any(0)
    if not present.any():
        raise ValueError("Batch without labels")
    return (per_head[present] * weights[present]).sum() / weights[present].sum()


def interval_multitask_loss(prediction, truth, weights, lower, upper, denominator, aux_kind):
    """Official soft-threshold interval error for pIC50; explicit point loss for auxiliaries."""
    mask = torch.isfinite(truth)
    potency = prediction[:, :4]
    distance = torch.relu(torch.nan_to_num(lower) - potency)
    distance = distance + torch.relu(potency - torch.nan_to_num(upper))
    distance = distance / denominator
    auxiliary = prediction[:, 4:] - torch.nan_to_num(truth[:, 4:])
    auxiliary = auxiliary.square() if aux_kind == "mse" else auxiliary.abs()
    error = torch.cat([distance, auxiliary], dim=1)
    per_head = (error * mask).sum(0) / mask.sum(0).clamp_min(1)
    present = mask.any(0)
    assert present.any()
    return (per_head[present] * weights[present]).sum() / weights[present].sum()


def inner_pic50_score(scorer, prediction, truth, lower, upper):
    scores = []
    for head in range(4):
        keep = np.isfinite(truth[:, head])
        assert keep.any()
        assert np.isfinite(prediction[keep, head]).all()
        assert np.isfinite(lower[keep, head]).all() and np.isfinite(upper[keep, head]).all()
        scores.append(
            float(
                scorer.rae_soft_threshold_absolute_error(
                    truth[keep, head],
                    prediction[keep, head],
                    y_true_lower=lower[keep, head],
                    y_true_upper=upper[keep, head],
                )
            )
        )
    assert np.isfinite(scores).all()
    return float(np.mean(scores)), scores


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = dict(
        seed=42,
        encoder_lr=1e-5,
        head_lr=1e-4,
        weight_decay=1e-4,
        max_epochs=60,
        patience=10,
        batch_size=64,
        inference_batch_size=128,
        dropout=0.1,
        inner_fraction=0.1,
        inner_group_column="group",
        device="cuda",
        num_threads=4,
        pIC50_weight=1.0,
        log2fc_weight=1.0,
        checkpoint_file="best.pt",
        training_loss="mse",
        auxiliary_loss="mse",
        checkpoint_selection="multitask_mse",
    )
    cfg.update(json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text()))
    assert cfg["pIC50_weight"] > 0 and cfg["log2fc_weight"] > 0
    assert cfg["training_loss"] in {"mse", "mae", "pic50_interval"}
    assert cfg["auxiliary_loss"] in {"mse", "mae"}
    assert cfg["checkpoint_selection"] in {"multitask_mse", "pic50_st_rae"}
    (out / "effective-config.json").write_text(json.dumps(cfg, indent=2))
    compounds = pl.read_parquet(input_path(inp, "compounds.parquet")).sort("compound_id")
    assignments = pl.read_parquet(input_path(inp, "parent-assignments.parquet")).sort("compound_id")
    labels = pl.read_parquet(input_path(inp, "measurements.parquet")).filter(
        pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition")
    )
    auxiliary = pl.read_parquet(
        inp / cfg["auxiliary_file"]
        if cfg.get("auxiliary_file")
        else input_path(inp, "targets.parquet")
    )
    endpoints = sorted(labels["endpoint"].unique())
    aux_heads = sorted(auxiliary["head"].unique())
    assert len(endpoints) == len(aux_heads) == 4
    assert {ep.split("_")[0] for ep in endpoints} == {head.split("@")[0] for head in aux_heads}
    heads = endpoints + aux_heads
    extra = None
    if cfg.get("extra_targets_file"):
        extra = pl.read_parquet(inp / cfg["extra_targets_file"])
        extra_heads = sorted(extra["head"].unique())
        assert not set(extra_heads) & set(heads)
        assert extra_heads and cfg.get("extra_weight", 1.0) > 0
        heads += extra_heads
    targets = pl.concat(
        [
            labels.select("compound_id", pl.col("endpoint").alias("head"), "value"),
            auxiliary.select("compound_id", "head", "value"),
        ]
    )
    if extra is not None:
        targets = pl.concat([targets, extra.select("compound_id", "head", "value")])
    assert targets.select("compound_id", "head").unique().height == targets.height
    assert set(targets["compound_id"]).issubset(set(assignments["compound_id"]))
    assert assignments["compound_id"].n_unique() == assignments.height
    assert compounds["compound_id"].n_unique() == compounds.height
    h = assignments.join(
        targets.pivot(on="head", index="compound_id", values="value"),
        on="compound_id",
        how="left",
        validate="1:1",
    ).sort("compound_id")
    assert h[cfg["inner_group_column"]].null_count() == 0
    y = h.select(heads).to_numpy().astype(np.float32)
    scorer, bounds_low, bounds_high = None, None, None
    if cfg["checkpoint_selection"] == "pic50_st_rae" or cfg["training_loss"] == "pic50_interval":
        spec = importlib.util.spec_from_file_location(
            "official", input_path(inp, "custom_scoring_functions.py")
        )
        scorer = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(scorer)
        bound_tables = [
            h.select("compound_id").join(
                labels.filter(pl.col("endpoint") == endpoint).select(
                    "compound_id", "conf_low", "conf_high"
                ),
                on="compound_id",
                how="left",
                validate="1:1",
                maintain_order="left",
            )
            for endpoint in endpoints
        ]
        bounds_low = np.column_stack([table["conf_low"].to_numpy() for table in bound_tables])
        bounds_high = np.column_stack([table["conf_high"].to_numpy() for table in bound_tables])
        measured = np.isfinite(y[:, :4])
        assert np.isfinite(bounds_low[measured]).all() and np.isfinite(bounds_high[measured]).all()
        assert np.all(bounds_low[measured] <= bounds_high[measured])
    labelled = np.isfinite(y).any(1)
    id_to_graph = {key: i for i, key in enumerate(compounds["compound_id"])}
    graph_indices = np.array([id_to_graph[key] for key in h["compound_id"]])
    # Trusted, explicitly declared own-run checkpoint; only MP state is transferred.
    checkpoint = torch.load(
        input_path(inp, cfg["checkpoint_file"]), map_location="cpu", weights_only=False
    )
    encoder_state = {
        key.removeprefix("mp."): value
        for key, value in checkpoint["state_dict"].items()
        if key.startswith("mp.")
    }
    assert encoder_state
    architecture = encoder_architecture(checkpoint)
    (out / "encoder-architecture.json").write_text(json.dumps(architecture, indent=2))
    torch.set_num_threads(cfg["num_threads"])
    torch.backends.cuda.matmul.allow_tf32 = False
    preparation_only = cfg.get("preparation_only", False)
    device = "cpu" if preparation_only else cfg["device"]
    weights = torch.tensor(
        [
            float(cfg.get("pic50_endpoint_weights", {}).get(ep.split("_")[0], cfg["pIC50_weight"]))
            for ep in endpoints
        ]
        + [cfg["log2fc_weight"]] * 4
        + [cfg.get("extra_weight", 1.0)] * (len(heads) - 8),
        device=device,
        dtype=torch.float32,
    )
    f = featurizers.SimpleMoleculeMolGraphFeaturizer()
    graphs = []
    for smiles in [] if preparation_only else compounds["standardized_smiles"]:
        molecule = Chem.MolFromSmiles(smiles)
        assert molecule is not None
        graphs.append(f(molecule))

    def graph(indices):
        batch = BatchMolGraph([graphs[int(i)] for i in indices])
        batch.to(device)
        return batch

    rows, fold_audits = [], []
    folds, anchor = h["fold"].to_numpy(), h["is_anchor"].to_numpy()
    groups = h[cfg["inner_group_column"]].to_numpy()
    for fold in sorted(h["fold"].unique()):
        rng = np.random.default_rng(cfg["seed"] + fold)
        torch.manual_seed(cfg["seed"] + fold)
        outer = (folds == fold) & ~anchor
        eligible = ~outer & labelled
        anchor_groups = set(groups[eligible & anchor])
        available_groups = sorted(set(groups[eligible]) - anchor_groups)
        if len(available_groups) < 2:
            raise ValueError("Too few nonanchor groups for inner validation")
        shuffled = rng.permutation(available_groups)
        val_groups = set(shuffled[: max(1, int(len(shuffled) * cfg["inner_fraction"]))])
        inner_val = eligible & np.isin(groups, list(val_groups))
        inner_train = eligible & ~inner_val
        assert not np.any(outer & (inner_val | inner_train))
        assert not set(groups[inner_train]) & set(groups[inner_val])
        assert np.all(inner_train[anchor & labelled])
        training, validation = np.flatnonzero(inner_train), np.flatnonzero(inner_val)
        means, stds = np.nanmean(y[training], axis=0), np.nanstd(y[training], axis=0)
        assert np.isfinite(means).all() and np.isfinite(stds).all()
        stds = np.where(stds < 1e-6, 1.0, stds)
        z = (y - means) / stds
        assert np.isfinite(y[validation]).sum(0).min() > 0
        directory = out / f"fold{fold}"
        directory.mkdir()
        interval_context = None
        if cfg["training_loss"] == "pic50_interval":
            # Scales, confidence bounds, and denominator use inner training only.
            low_z = (bounds_low - means[:4]) / stds[:4]
            high_z = (bounds_high - means[:4]) / stds[:4]
            denominators, verification = [], []
            for j in range(4):
                keep = np.isfinite(y[training, j])
                raw_y = y[training, j][keep]
                lo, hi = bounds_low[training, j][keep], bounds_high[training, j][keep]
                mean = np.mean(raw_y)
                denominator = np.mean(np.maximum(lo - mean, 0) + np.maximum(mean - hi, 0))
                assert denominator > 0 and np.isfinite(denominator)
                denominators.append(denominator / stds[j])
                probe = np.full(len(raw_y), float(mean) + 0.123 * float(stds[j]))
                ratio = np.mean(np.maximum(lo - probe, 0) + np.maximum(probe - hi, 0)) / denominator
                official = scorer.rae_soft_threshold_absolute_error(
                    raw_y, probe, y_true_lower=lo, y_true_upper=hi
                )
                assert np.isclose(ratio, official, atol=1e-12)
                verification.append(
                    dict(
                        head=heads[j],
                        training_n=len(raw_y),
                        raw_denominator=denominator,
                        equivalent_official_probe_score=float(official),
                    )
                )
            interval_context = (
                low_z,
                high_z,
                torch.tensor(denominators, dtype=torch.float32, device=device),
            )
            (directory / "interval-training-context.json").write_text(
                json.dumps(
                    dict(
                        normalization="Per-head mean soft error of inner-training mean predictor",
                        auxiliary_loss=cfg["auxiliary_loss"],
                        verification=verification,
                    ),
                    indent=2,
                )
            )
            h.select("compound_id").filter(pl.Series(inner_train)).with_columns(
                *[pl.Series(heads[j] + "_conf_low", bounds_low[training, j]) for j in range(4)],
                *[pl.Series(heads[j] + "_conf_high", bounds_high[training, j]) for j in range(4)],
            ).write_parquet(directory / "inner-training-confidence-bounds.parquet")
        role = np.where(
            outer,
            "outer_holdout",
            np.where(inner_val, "inner_validation", np.where(inner_train, "training", "no_labels")),
        )
        h.select("compound_id", "fold", "is_anchor", cfg["inner_group_column"]).with_columns(
            pl.Series("role", role),
            pl.lit(fold).alias("outer_fold"),
            pl.lit(cfg["split_run"]).alias("split_id"),
        ).write_parquet(directory / "membership.parquet")
        if preparation_only:
            fold_audits.append(
                dict(
                    fold=fold,
                    training=len(training),
                    inner_validation=len(validation),
                    outer_compounds=int(outer.sum()),
                    training_labels=dict(
                        zip(heads, np.isfinite(y[training]).sum(0).tolist(), strict=True)
                    ),
                    validation_labels=dict(
                        zip(heads, np.isfinite(y[validation]).sum(0).tolist(), strict=True)
                    ),
                    outer_labels_excluded=dict(
                        zip(heads, np.isfinite(y[outer]).sum(0).tolist(), strict=True)
                    ),
                    means=means.tolist(),
                    stds=stds.tolist(),
                )
            )
            continue
        model = MultitaskModel(encoder_state, cfg["dropout"], len(heads), architecture).to(device)
        optimizer = torch.optim.AdamW(
            [
                {"params": model.mp.parameters(), "lr": cfg["encoder_lr"]},
                {"params": model.head.parameters(), "lr": cfg["head_lr"]},
            ],
            weight_decay=cfg["weight_decay"],
        )
        best, best_epoch, stale = float("inf"), -1, 0
        history = []
        for epoch in range(cfg["max_epochs"]):
            model.train()
            total, count = 0.0, 0
            order = rng.permutation(training)
            for start in range(0, len(order), cfg["batch_size"]):
                batch = order[start : start + cfg["batch_size"]]
                optimizer.zero_grad(set_to_none=True)
                prediction = model(graph(graph_indices[batch]))
                truth = torch.as_tensor(z[batch], device=device)
                if interval_context is None:
                    loss = masked_loss(prediction, truth, weights, cfg["training_loss"])
                else:
                    lower, upper, denominator = interval_context
                    loss = interval_multitask_loss(
                        prediction,
                        truth,
                        weights,
                        torch.as_tensor(lower[batch], dtype=torch.float32, device=device),
                        torch.as_tensor(upper[batch], dtype=torch.float32, device=device),
                        denominator,
                        cfg["auxiliary_loss"],
                    )
                if not torch.isfinite(loss):
                    raise ValueError("Nonfinite training loss")
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0, error_if_nonfinite=True)
                optimizer.step()
                total += float(loss.detach()) * len(batch)
                count += len(batch)
            model.eval()
            with torch.inference_mode():
                pv = torch.cat(
                    [
                        model(
                            graph(
                                graph_indices[
                                    validation[start : start + cfg["inference_batch_size"]]
                                ]
                            )
                        )
                        for start in range(0, len(validation), cfg["inference_batch_size"])
                    ]
                )
                val_loss = float(
                    masked_loss(pv, torch.as_tensor(z[validation], device=device), weights)
                )
                per_head = np.nanmean((pv.cpu().numpy() - z[validation]) ** 2, axis=0)
            assert np.isfinite(val_loss)
            selection_score = val_loss
            record = dict(
                fold=fold,
                epoch=epoch,
                validation_mse=val_loss,
                validation_head_mse=dict(zip(heads, per_head.astype(float), strict=True)),
            )
            record["train_" + cfg["training_loss"]] = total / count
            if cfg["checkpoint_selection"] == "pic50_st_rae":
                selection_score, st_by_head = inner_pic50_score(
                    scorer,
                    (pv.cpu().numpy() * stds + means)[:, :4],
                    y[validation, :4],
                    bounds_low[validation],
                    bounds_high[validation],
                )
                record["validation_pic50_st_rae"] = selection_score
                selection_indices = [
                    j
                    for j, ep in enumerate(endpoints)
                    if ep.split("_")[0]
                    in cfg.get("selection_targets", [ep.split("_")[0] for ep in endpoints])
                ]
                assert selection_indices
                selection_score = float(np.mean([st_by_head[j] for j in selection_indices]))
                record["validation_focus_st_rae"] = selection_score
                record["validation_pic50_head_st_rae"] = dict(
                    zip(endpoints, st_by_head, strict=True)
                )
            record["selection_score"] = selection_score
            history.append(record)
            with (directory / "history.jsonl").open("a") as stream:
                stream.write(json.dumps(record) + "\n")
            print(json.dumps(record), flush=True)
            state = dict(
                state_dict=model.state_dict(),
                optimizer=optimizer.state_dict(),
                epoch=epoch,
                heads=heads,
                means=means.tolist(),
                stds=stds.tolist(),
                config=cfg,
                encoder_architecture=architecture,
                fold=fold,
                split_id=cfg["split_run"],
            )
            torch.save(state, directory / "last.pt")
            if selection_score < best:
                best, best_epoch, stale = selection_score, epoch, 0
                best_state = copy.deepcopy(model.state_dict())
                torch.save(state, directory / "best.pt")
            else:
                stale += 1
            if stale >= cfg["patience"]:
                break
        model.load_state_dict(best_state, strict=True)
        model.eval()
        embeddings, predictions = [], []
        with torch.inference_mode():
            for start in range(0, compounds.height, cfg["inference_batch_size"]):
                batch = np.arange(start, min(start + cfg["inference_batch_size"], compounds.height))
                emb = model.embed(graph(batch))
                embeddings.append(emb.cpu().numpy())
                predictions.append(model.head(emb).cpu().numpy() * stds + means)
        embeddings, predictions = np.concatenate(embeddings), np.concatenate(predictions)
        assert embeddings.shape == (compounds.height, 600) and np.isfinite(embeddings).all()
        assert np.isfinite(predictions).all()
        pl.DataFrame(embeddings, schema=[f"emb_{i:03d}" for i in range(600)]).insert_column(
            0, compounds["compound_id"]
        ).with_columns(
            pl.lit(fold).alias("fold"), pl.lit(cfg["split_run"]).alias("split_id")
        ).write_parquet(directory / "external-embeddings.parquet")
        pl.DataFrame(predictions, schema=heads).insert_column(
            0, compounds["compound_id"]
        ).with_columns(
            pl.lit(fold).alias("fold"), pl.lit(cfg["split_run"]).alias("split_id")
        ).write_parquet(directory / "direct-predictions.parquet")
        outer_ids = h.filter(pl.Series(outer))["compound_id"]
        for j, endpoint in enumerate(endpoints):
            evaluated = labels.filter(
                (pl.col("endpoint") == endpoint) & pl.col("compound_id").is_in(outer_ids.implode())
            )
            evaluated = evaluated.join(assignments, on="compound_id", validate="1:1").sort(
                "compound_id"
            )
            prediction = predictions[[id_to_graph[key] for key in evaluated["compound_id"]], j]
            rows.append(
                evaluated.select(
                    "compound_id",
                    "endpoint",
                    "fold",
                    "is_neighbor",
                    "conf_low",
                    "conf_high",
                    pl.col("value").alias("y_true"),
                ).with_columns(
                    pl.Series("prediction", prediction.astype(float)),
                    pl.lit(cfg.get("method", "external_multitask_ft_direct")).alias("method"),
                    pl.lit(cfg["split_run"]).alias("split_id"),
                )
            )
        fold_audits.append(
            dict(
                fold=fold,
                training=len(training),
                inner_validation=len(validation),
                outer_compounds=int(outer.sum()),
                best_epoch=best_epoch,
                best_validation_mse=history[best_epoch]["validation_mse"],
                best_selection_score=best,
                checkpoint_selection=cfg["checkpoint_selection"],
                epochs=epoch + 1,
                training_labels=dict(
                    zip(heads, np.isfinite(y[training]).sum(0).tolist(), strict=True)
                ),
                validation_labels=dict(
                    zip(heads, np.isfinite(y[validation]).sum(0).tolist(), strict=True)
                ),
                outer_labels_excluded=dict(
                    zip(heads, np.isfinite(y[outer]).sum(0).tolist(), strict=True)
                ),
            )
        )
    if preparation_only:
        (out / "summary.json").write_text(
            json.dumps(
                dict(
                    config=cfg,
                    heads=heads,
                    folds=fold_audits,
                    preparation_only=True,
                    training_executed=False,
                    external_encoder_tensors=len(encoder_state),
                    encoder_architecture=architecture,
                ),
                indent=2,
            )
        )
        return
    oof = pl.concat(rows)
    expected = labels.join(assignments, on="compound_id", validate="m:1").filter(
        ~pl.col("is_anchor")
    )
    assert oof.height == expected.height
    assert oof.select("compound_id", "endpoint").unique().height == oof.height
    if "expected_oof_rows" in cfg:
        assert oof.height == cfg["expected_oof_rows"]
    oof.write_parquet(out / "finetune-oof.parquet")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                config=cfg,
                heads=heads,
                folds=fold_audits,
                oof_rows=oof.height,
                dimensions=600,
                challenge_rows=compounds.height,
                selection="group-held-out inner loss; all outer compound labels excluded",
                external_checkpoint_epoch=checkpoint.get("epoch"),
                encoder_architecture=architecture,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
