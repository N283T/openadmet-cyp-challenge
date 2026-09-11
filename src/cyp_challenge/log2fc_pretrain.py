"""Fold-safe auxiliary pretraining from CheMeleon, with explicit best-checkpoint selection."""

import copy
import json
import math
import os
import time
from pathlib import Path

import numpy as np
import polars as pl
import torch
from chemprop import featurizers
from chemprop import nn as cnn
from chemprop.data import BatchMolGraph
from rdkit import Chem
from scipy.stats import spearmanr
from torch import nn

from .storage import digest, write_json


def partition(ids, folds, outer_fold, seed, fraction, excluded_folds=None):
    excluded = [outer_fold] if excluded_folds is None else excluded_folds
    if outer_fold not in excluded:
        raise ValueError("The outer fold must be excluded")
    eligible = np.flatnonzero(~np.isin(folds, excluded))
    heldout = np.flatnonzero(np.isin(folds, excluded))
    shuffled = np.random.default_rng(seed).permutation(eligible)
    n_val = max(1, int(len(eligible) * fraction))
    validation, training = shuffled[:n_val], shuffled[n_val:]
    if not len(training) or not len(heldout):
        raise ValueError("Empty training or outer held-out set")
    assert not (set(ids[i] for i in training) & set(ids[i] for i in heldout))
    return training, validation, heldout


def target_scaler(y, training):
    means = np.nanmean(y[training], axis=0).astype(np.float32)
    stds = np.nanstd(y[training], axis=0).astype(np.float32)
    if not np.isfinite(means).all() or not np.isfinite(stds).all():
        raise ValueError("A target has no training labels")
    return means, np.where(stds < 1e-6, 1, stds)


class Network(nn.Module):
    def __init__(self, checkpoint, config, n_heads):
        super().__init__()
        self.mp = cnn.BondMessagePassing(**checkpoint["hyper_parameters"])
        if config["initialization"] == "chemeleon":
            self.mp.load_state_dict(checkpoint["state_dict"], strict=True)
        elif config["initialization"] != "scratch":
            raise ValueError("Unknown initialization")
        self.agg = cnn.MeanAggregation()
        self.bn = nn.BatchNorm1d(self.mp.output_dim)
        self.head = nn.Sequential(
            nn.Linear(self.mp.output_dim, config["ffn_hidden_dim"]),
            nn.ReLU(),
            nn.Dropout(config["ffn_dropout"]),
            nn.Linear(config["ffn_hidden_dim"], n_heads),
        )

    def fingerprint(self, graph):
        return self.bn(self.agg(self.mp(graph), graph.batch))

    def forward(self, graph):
        return self.head(self.fingerprint(graph))


def main():
    inputs, output = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
    seed = config["seed"]
    torch.manual_seed(seed)
    torch.set_num_threads(4)
    torch.use_deterministic_algorithms(True)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    device = torch.device("cuda")
    checkpoint_path = next(inputs.glob("*-chemeleon_mp.pt"))
    assert digest(checkpoint_path) == config["checkpoint_sha256"]
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
    provenance = json.loads(next(inputs.glob("*-provenance.json")).read_text())
    targets = pl.read_parquet(next(inputs.glob("*-targets.parquet")))
    compounds = pl.read_parquet(next(inputs.glob("*-compounds.parquet"))).sort("compound_id")
    assignments = pl.read_parquet(next(inputs.glob("*-assignments.parquet"))).filter(
        pl.col("strategy") == config.get("outer_strategy", "umap")
    )
    heads = sorted(targets["head"].unique())
    wide = (
        targets.pivot(on="head", index="compound_id", values="value")
        .join(
            assignments.select("compound_id", "fold"), on="compound_id", how="left", validate="1:1"
        )
        .join(compounds, on="compound_id", how="left", validate="1:1")
        .sort("compound_id")
    )
    assert wide["fold"].null_count() == wide["standardized_smiles"].null_count() == 0
    ids = wide["compound_id"].to_list()
    y = wide.select(heads).to_numpy().astype(np.float32)
    folds = wide["fold"].to_numpy()
    outer = config["outer_fold"]
    excluded = config.get("excluded_folds", [outer])
    train, validation, heldout = partition(
        ids, folds, outer, seed, config["inner_validation_fraction"], excluded_folds=excluded
    )
    means, stds = target_scaler(y, train)
    z = (y - means) / stds
    roles = ["outer_holdout"] * len(ids)
    for indexes, role in [(train, "inner_train"), (validation, "inner_validation")]:
        for i in indexes:
            roles[i] = role
    pl.DataFrame({"compound_id": ids, "outer_fold": outer, "role": roles}).write_parquet(
        output / "membership.parquet"
    )
    metadata = {
        "config": config,
        "heads": heads,
        "target_means": means.tolist(),
        "target_stds": stds.tolist(),
        "counts": {
            "train": len(train),
            "validation": len(validation),
            "outer_holdout": len(heldout),
        },
        "checkpoint_provenance": provenance,
        "architecture": checkpoint["hyper_parameters"],
        "torch": torch.__version__,
        "device": torch.cuda.get_device_name(),
        "cuda": torch.version.cuda,
        "embedding": "2048d mean pooled message passing followed by trained BatchNorm",
        "inner_validation": "seeded 90/10 within outer-training molecules, for stopping only",
    }
    write_json(output / "metadata.json", metadata)
    featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
    graph_by_id = {
        cid: featurizer(Chem.MolFromSmiles(smi))
        for cid, smi in compounds.select("compound_id", "standardized_smiles").iter_rows()
    }

    def graph(batch_ids):
        result = BatchMolGraph([graph_by_id[cid] for cid in batch_ids])
        result.to(device)
        return result

    model = Network(checkpoint, config, len(heads)).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=config["learning_rate"])
    batch_size = config["batch_size"]
    steps = math.ceil(len(train) / batch_size)
    warmup = config["warmup_epochs"] * steps
    total_steps = config["max_epochs"] * steps
    peak = config["learning_rate"] * config["lr_ratio"]
    final_lr = config["learning_rate"] * 0.1
    history, best, best_epoch, stale, step = [], float("inf"), -1, 0, 0
    rng = np.random.default_rng(seed)

    def predict(indices):
        model.eval()
        predictions = []
        with torch.inference_mode():
            for start in range(0, len(indices), batch_size):
                batch = indices[start : start + batch_size]
                predictions.append(model(graph([ids[i] for i in batch])).cpu().numpy())
        return np.concatenate(predictions)

    for epoch in range(config["max_epochs"]):
        started = time.monotonic()
        model.train()
        shuffled = rng.permutation(train)
        train_sum, train_count = 0.0, 0
        for batch in np.array_split(shuffled, steps):
            # array_split prevents a singleton BatchNorm training batch.
            lr = (
                config["learning_rate"] + (peak - config["learning_rate"]) * step / max(1, warmup)
                if step < warmup
                else peak * (final_lr / peak) ** ((step - warmup) / max(1, total_steps - warmup))
            )
            for group in optimizer.param_groups:
                group["lr"] = lr
            truth = torch.tensor(z[batch], device=device)
            mask = torch.isfinite(truth)
            optimizer.zero_grad(set_to_none=True)
            prediction = model(graph([ids[i] for i in batch]))
            loss = ((prediction[mask] - truth[mask]) ** 2).mean()
            if not torch.isfinite(loss):
                raise ValueError("Nonfinite training loss")
            loss.backward()
            optimizer.step()
            train_sum += float(loss.detach()) * int(mask.sum())
            train_count += int(mask.sum())
            step += 1
        p = predict(validation)
        val_loss = float(np.nanmean((p - z[validation]) ** 2))
        if not np.isfinite(val_loss):
            raise ValueError("Nonfinite validation loss")
        improved = val_loss < best
        if improved:
            best, best_epoch, stale = val_loss, epoch, 0
            state = {
                "state_dict": {k: v.detach().cpu() for k, v in model.state_dict().items()},
                "metadata": metadata,
                "epoch": epoch,
                "validation_mse_z": best,
            }
            torch.save(state, output / "best.pt")
        else:
            stale += 1
        row = {
            "epoch": epoch,
            "train_mse_z": train_sum / train_count,
            "val_mse_z": val_loss,
            "learning_rate": lr,
            "seconds": time.monotonic() - started,
            "best_epoch": best_epoch,
        }
        history.append(row)
        with (output / "history.jsonl").open("a") as stream:
            stream.write(json.dumps(row) + "\n")
        torch.save(
            {
                "state_dict": model.state_dict(),
                "optimizer": optimizer.state_dict(),
                "epoch": epoch,
                "step": step,
                "rng": copy.deepcopy(rng.bit_generator.state),
                "torch_rng": torch.get_rng_state(),
                "cuda_rng": torch.cuda.get_rng_state_all(),
                "metadata": metadata,
            },
            output / "last.pt",
        )
        print(json.dumps(row), flush=True)
        if stale >= config["patience"]:
            break
    best_state = torch.load(output / "best.pt", map_location="cpu", weights_only=False)
    model.load_state_dict(best_state["state_dict"], strict=True)
    outer_pred = predict(heldout) * stds + means
    rows, metrics = [], []
    for j, head in enumerate(heads):
        mask = np.isfinite(y[heldout, j])
        error = outer_pred[mask, j] - y[heldout[mask], j]
        metrics.append(
            {
                "head": head,
                "n": int(mask.sum()),
                "rmse": float(np.sqrt((error**2).mean())),
                "mae": float(np.abs(error).mean()),
                "spearman": float(spearmanr(y[heldout[mask], j], outer_pred[mask, j]).statistic),
            }
        )
        for k, i in enumerate(heldout):
            rows.append(
                {
                    "compound_id": ids[i],
                    "endpoint": head,
                    "split_id": config["split_run"],
                    "fold": outer,
                    "y_true": float(y[i, j]),
                    "prediction": float(outer_pred[k, j]),
                }
            )
    pl.DataFrame(rows).write_parquet(output / "outer_predictions.parquet")
    # Persist the selected encoder's features for later readout experiments; no readout fit here.
    model.eval()
    for start in range(0, compounds.height, batch_size):
        batch = compounds.slice(start, batch_size)
        with torch.inference_mode():
            g = graph(batch["compound_id"].to_list())
            if config.get("export_predictions_only", False):
                values = model(g).cpu().numpy() * stds + means
            else:
                values = model.fingerprint(g).cpu().numpy()
        assert np.isfinite(values).all()
        if config.get("export_predictions_only", False):
            pl.DataFrame(
                {
                    "compound_id": batch["compound_id"],
                    **{h: values[:, j] for j, h in enumerate(heads)},
                }
            ).write_parquet(output / f"predictions-{start:05d}.parquet")
        else:
            pl.DataFrame(
                {"compound_id": batch["compound_id"], "fold": outer, "embedding": values.tolist()}
            ).write_parquet(output / f"embeddings-{start:05d}.parquet")
    write_json(
        output / "summary.json",
        {
            **metadata,
            "best_epoch": best_epoch,
            "best_validation_mse_z": best,
            "epochs": len(history),
            "outer_metrics": metrics,
            "best_checkpoint_sha256": digest(output / "best.pt"),
            "feature_rows": compounds.height,
            "dimension": model.mp.output_dim,
        },
    )
    print("COMPLETE", json.dumps(metrics), flush=True)


if __name__ == "__main__":
    main()
