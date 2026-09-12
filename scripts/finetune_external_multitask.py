"""Fold-safe external600 -> four pIC50 plus four actual log2FC challenge heads.

Outer nonanchor compounds lose ALL eight labels during fitting and checkpoint
selection. Anchors stay in training. Inner validation separates supplied molecular
groups; groups containing anchors cannot enter inner validation.
"""

import copy
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


class MultitaskModel(torch.nn.Module):
    def __init__(self, encoder_state, dropout):
        super().__init__()
        self.mp = cnn.BondMessagePassing(d_h=600, depth=5)
        self.mp.load_state_dict(encoder_state, strict=True)
        self.aggregation = cnn.MeanAggregation()
        self.head = torch.nn.Sequential(
            torch.nn.Linear(600, 256),
            torch.nn.ReLU(),
            torch.nn.Dropout(dropout),
            torch.nn.Linear(256, 8),
        )

    def embed(self, graph):
        return self.aggregation(self.mp(graph), graph.batch)

    def forward(self, graph):
        return self.head(self.embed(graph))


def masked_loss(prediction, truth, weights):
    mask = torch.isfinite(truth)
    per_head = ((prediction - torch.nan_to_num(truth)).square() * mask).sum(0)
    per_head = per_head / mask.sum(0).clamp_min(1)
    present = mask.any(0)
    if not present.any():
        raise ValueError("Batch without labels")
    return (per_head[present] * weights[present]).sum() / weights[present].sum()


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
    )
    cfg.update(json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text()))
    assert cfg["pIC50_weight"] > 0 and cfg["log2fc_weight"] > 0
    (out / "effective-config.json").write_text(json.dumps(cfg, indent=2))
    compounds = pl.read_parquet(input_path(inp, "compounds.parquet")).sort("compound_id")
    assignments = pl.read_parquet(input_path(inp, "parent-assignments.parquet")).sort("compound_id")
    labels = pl.read_parquet(input_path(inp, "measurements.parquet")).filter(
        pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition")
    )
    auxiliary = pl.read_parquet(input_path(inp, "targets.parquet"))
    endpoints = sorted(labels["endpoint"].unique())
    aux_heads = sorted(auxiliary["head"].unique())
    assert len(endpoints) == len(aux_heads) == 4
    assert {ep.split("_")[0] for ep in endpoints} == {head.split("@")[0] for head in aux_heads}
    heads = endpoints + aux_heads
    targets = pl.concat(
        [
            labels.select("compound_id", pl.col("endpoint").alias("head"), "value"),
            auxiliary.select("compound_id", "head", "value"),
        ]
    )
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
    torch.set_num_threads(cfg["num_threads"])
    torch.backends.cuda.matmul.allow_tf32 = False
    preparation_only = cfg.get("preparation_only", False)
    device = "cpu" if preparation_only else cfg["device"]
    weights = torch.tensor(
        [cfg["pIC50_weight"]] * 4 + [cfg["log2fc_weight"]] * 4, device=device, dtype=torch.float32
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
        model = MultitaskModel(encoder_state, cfg["dropout"]).to(device)
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
                loss = masked_loss(
                    model(graph(graph_indices[batch])),
                    torch.as_tensor(z[batch], device=device),
                    weights,
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
            record = dict(
                fold=fold,
                epoch=epoch,
                train_mse=total / count,
                validation_mse=val_loss,
                validation_head_mse=dict(zip(heads, per_head.astype(float), strict=True)),
            )
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
                fold=fold,
                split_id=cfg["split_run"],
            )
            torch.save(state, directory / "last.pt")
            if val_loss < best:
                best, best_epoch, stale = val_loss, epoch, 0
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
                best_validation_mse=best,
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
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
