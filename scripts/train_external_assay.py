"""Assay-separated scratch Chemprop600 external pretraining; no challenge labels.

Continuous interval regression and ranking share assay potency outputs; binary
classification uses separate logits. Select objective then assay uniformly.
Ranking pairs never cross assays and require nonoverlapping censor intervals.
"""

import copy
import json
import math
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
import torch.nn.functional as F
from chemprop import featurizers
from chemprop import nn as cnn
from chemprop.data import BatchMolGraph
from external_loss_utils import interval_mse, ordered_pairs, ranking_loss
from rdkit import Chem


def read_input(root, name):
    paths = [p for p in root.iterdir() if p.name.partition("-")[2] == name]
    if len(paths) != 1:
        raise ValueError(f"Expected one input {name}; got {paths}")
    return pl.read_parquet(paths[0])


class AssayModel(torch.nn.Module):
    def __init__(self, n_assays, dropout):
        super().__init__()
        self.mp = cnn.BondMessagePassing(d_h=600, depth=5)
        self.aggregation = cnn.MeanAggregation()
        self.head = torch.nn.Sequential(
            torch.nn.Linear(600, 600),
            torch.nn.ReLU(),
            torch.nn.Dropout(dropout),
            torch.nn.Linear(600, n_assays * 2),
        )

    def embed(self, graph):
        return self.aggregation(self.mp(graph), graph.batch)

    def forward(self, graph):
        return self.head(self.embed(graph))


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = dict(
        seed=42,
        max_epochs=60,
        patience=10,
        batch_size=128,
        dropout=0.1,
        learning_rate=1e-3,
        min_learning_rate=1e-4,
        weight_decay=1e-4,
        validation_fraction=0.1,
        min_train_observations=8,
        rank_min_gap=0.25,
        max_train_pairs=4096,
        max_validation_pairs=2048,
        steps_per_epoch=200,
        device="cuda",
        num_threads=4,
        interval_regression=True,
        continuous_weight=1.0,
        ranking_weight=1.0,
        classification_weight=1.0,
    )
    cfg.update(json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text()))
    modes = {
        "continuous": ["continuous"],
        "ranking": ["ranking"],
        "classification": ["classification"],
        "mixed": ["continuous", "ranking", "classification"],
    }
    objectives = modes[cfg["mode"]]
    (out / "effective-config.json").write_text(json.dumps(cfg, indent=2))
    compounds = read_input(inp, "external-compounds.parquet").sort("compound_id")
    observations = read_input(inp, "observations.parquet")
    challenge = read_input(inp, "compounds.parquet").sort("compound_id")
    assert compounds["compound_id"].n_unique() == compounds.height
    assert challenge["compound_id"].n_unique() == challenge.height
    assert not set(compounds["compound_id"]) & set(challenge["compound_id"])
    if "connectivity" in challenge.columns:
        assert not set(compounds["connectivity"]) & set(challenge["connectivity"])
    # Exact replicates use median potency. Mixed/censored intervals intersect;
    # inconsistent intervals and contradictory binary replicates are masked.
    if observations.select("compound_id", "assay_id").unique().height != observations.height:
        observations = observations.group_by("compound_id", "assay_id").agg(
            pl.col("value").median(),
            pl.when(pl.col("binary_label").drop_nulls().n_unique() == 1)
            .then(pl.col("binary_label").drop_nulls().first())
            .otherwise(None)
            .alias("binary_label"),
            pl.col("source").first(),
            pl.when((pl.col("lower_bound") == pl.col("upper_bound")).fill_null(False).all())
            .then(pl.col("lower_bound").median())
            .otherwise(pl.col("lower_bound").max())
            .alias("lower_bound"),
            pl.when((pl.col("lower_bound") == pl.col("upper_bound")).fill_null(False).all())
            .then(pl.col("upper_bound").median())
            .otherwise(pl.col("upper_bound").min())
            .alias("upper_bound"),
        )
        observations = observations.with_columns(
            pl.when(pl.col("lower_bound") > pl.col("upper_bound"))
            .then(None)
            .otherwise(pl.col(name))
            .alias(name)
            for name in ["lower_bound", "upper_bound"]
        )
    assert set(observations["compound_id"]).issubset(set(compounds["compound_id"]))
    indices = {key: i for i, key in enumerate(compounds["compound_id"])}
    rng = np.random.default_rng(cfg["seed"])
    val_rng = np.random.default_rng(cfg["seed"] + 739)
    keys = rng.permutation(sorted(compounds["connectivity"].unique()))
    val_keys = set(keys[: max(1, int(len(keys) * cfg["validation_fraction"]))])
    is_val = np.array([key in val_keys for key in compounds["connectivity"]])
    compounds.select("compound_id", "connectivity").with_columns(
        pl.Series("role", np.where(is_val, "validation", "training"))
    ).write_parquet(out / "membership.parquet")
    tasks, task_meta = [], []
    assays = sorted(observations["assay_id"].unique())
    for assay_index, assay in enumerate(assays):
        rows = observations.filter(pl.col("assay_id") == assay).sort("compound_id")
        ids = np.array([indices[key] for key in rows["compound_id"]])
        value = rows["value"].cast(pl.Float64).to_numpy()
        binary = rows["binary_label"].cast(pl.Float64).to_numpy()
        assert np.isin(binary[np.isfinite(binary)], [0, 1]).all()
        if "lower_bound" in rows.columns:
            lower = rows["lower_bound"].cast(pl.Float64).to_numpy()
            upper = rows["upper_bound"].cast(pl.Float64).to_numpy()
            known = np.isfinite(lower) | np.isfinite(upper)
            lower = np.where(known & np.isnan(lower), -np.inf, lower)
            upper = np.where(known & np.isnan(upper), np.inf, upper)
        else:
            rel = rows["relation"].to_numpy() if "relation" in rows.columns else "="
            lower = np.where(np.isin(rel, ["=", "==", ">", ">="]), value, -np.inf)
            upper = np.where(np.isin(rel, ["=", "==", "<", "<="]), value, np.inf)
            lower = np.where(np.isfinite(value), lower, np.nan)
            upper = np.where(np.isfinite(value), upper, np.nan)
        exact = np.isfinite(lower) & (lower == upper)
        train_exact = exact & ~is_val[ids]
        mean = float(np.mean(lower[train_exact])) if train_exact.any() else 0.0
        std = float(np.std(lower[train_exact])) if train_exact.sum() > 1 else 1.0
        std = max(std, 0.1)
        known = (np.isfinite(lower) | np.isfinite(upper)) & (lower <= upper)
        assert not np.any(~np.isnan(lower) & ~np.isnan(upper) & (lower > upper))
        for objective in objectives:
            eligible = np.isfinite(binary) if objective == "classification" else known
            if objective == "continuous":
                if not cfg["interval_regression"]:
                    eligible = exact
                if train_exact.sum() < 2:
                    continue
            tr = np.flatnonzero(eligible & ~is_val[ids])
            va = np.flatnonzero(eligible & is_val[ids])
            if len(tr) < cfg["min_train_observations"]:
                continue
            if objective == "classification" and len(np.unique(binary[tr])) != 2:
                continue
            tp = vp = None
            if objective == "ranking":
                tp = ordered_pairs(
                    lower[tr], upper[tr], cfg["max_train_pairs"], cfg["rank_min_gap"], rng
                )
                vp = ordered_pairs(
                    lower[va], upper[va], cfg["max_validation_pairs"], cfg["rank_min_gap"], val_rng
                )
                if not len(tp):
                    continue
            task = dict(
                assay_id=assay,
                source=rows["source"][0],
                objective=objective,
                column=assay_index * 2 + int(objective == "classification"),
                ids=ids,
                train=tr,
                validation=va,
                lower=(lower - mean) / std,
                upper=(upper - mean) / std,
                binary=binary,
                train_pairs=tp,
                validation_pairs=vp,
            )
            tasks.append(task)
            task_meta.append(
                dict(
                    assay_id=assay,
                    source=rows["source"][0],
                    objective=objective,
                    column=task["column"],
                    mean=mean,
                    std=std,
                    train_n=len(tr),
                    validation_n=len(va),
                    train_pairs=None if tp is None else len(tp),
                    validation_pairs=None if vp is None else len(vp),
                )
            )
    if not tasks:
        raise ValueError("No eligible training tasks")
    # Allocate outputs only for eligible assays, not every raw assay record.
    old_assays = assays
    assays = sorted({task["assay_id"] for task in tasks})
    assay_lookup = {assay: i for i, assay in enumerate(assays)}
    for task, meta in zip(tasks, task_meta, strict=True):
        task["column"] = assay_lookup[task["assay_id"]] * 2 + int(
            task["objective"] == "classification"
        )
        meta["column"] = task["column"]
    print(
        json.dumps(dict(raw_assays=len(old_assays), eligible_assays=len(assays), tasks=len(tasks))),
        flush=True,
    )
    (out / "tasks.json").write_text(json.dumps(task_meta, indent=2))
    torch.manual_seed(cfg["seed"])
    torch.set_num_threads(cfg["num_threads"])
    torch.backends.cuda.matmul.allow_tf32 = False
    device = cfg["device"]
    featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
    graphs = []
    for smiles in compounds["standardized_smiles"]:
        molecule = Chem.MolFromSmiles(smiles)
        assert molecule is not None
        graphs.append(featurizer(molecule))

    def graph(batch):
        result = BatchMolGraph([graphs[int(i)] for i in batch])
        result.to(device)
        return result

    model = AssayModel(len(assays), cfg["dropout"]).to(device)
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=cfg["learning_rate"], weight_decay=cfg["weight_decay"]
    )

    def objective_loss(task, positions, pairs=None, all_predictions=None):
        selected = positions[pairs.flatten()] if pairs is not None else positions
        if all_predictions is None:
            unique, inverse = np.unique(task["ids"][selected], return_inverse=True)
            prediction = model(graph(unique))[:, task["column"]]
            prediction = prediction[torch.as_tensor(inverse, device=device)]
        else:
            prediction = torch.as_tensor(
                all_predictions[np.searchsorted(val_ids, task["ids"][selected]), task["column"]],
                device=device,
            )
        if task["objective"] == "ranking":
            prediction = prediction.reshape(-1, 2)
            return ranking_loss(prediction[:, 0], prediction[:, 1])
        if task["objective"] == "classification":
            return F.binary_cross_entropy_with_logits(
                prediction,
                torch.as_tensor(task["binary"][selected], dtype=torch.float32, device=device),
            )
        return interval_mse(
            prediction,
            torch.as_tensor(task["lower"][selected], dtype=torch.float32, device=device),
            torch.as_tensor(task["upper"][selected], dtype=torch.float32, device=device),
        )

    val_ids = np.flatnonzero(is_val)
    available = {name: [t for t in tasks if t["objective"] == name] for name in objectives}
    available = {name: ts for name, ts in available.items() if ts}
    by_source = {
        name: {
            source: [t for t in ts if t["source"] == source]
            for source in sorted({t["source"] for t in ts})
        }
        for name, ts in available.items()
    }
    best, stale, best_epoch = float("inf"), 0, -1
    for epoch in range(cfg["max_epochs"]):
        phase = max(0, epoch - 3) / max(1, cfg["max_epochs"] - 4)
        lr = (
            cfg["learning_rate"] * (epoch + 1) / 3
            if epoch < 3
            else (
                cfg["min_learning_rate"]
                + (cfg["learning_rate"] - cfg["min_learning_rate"])
                * (1 + math.cos(math.pi * phase))
                / 2
            )
        )
        for group in optimizer.param_groups:
            group["lr"] = lr
        model.train()
        train_losses = {name: [] for name in available}
        for _ in range(cfg["steps_per_epoch"]):
            objective = rng.choice(list(available))
            source = rng.choice(list(by_source[objective]))
            pool_tasks = by_source[objective][source]
            task = pool_tasks[rng.integers(len(pool_tasks))]
            if objective == "ranking":
                pool = task["train_pairs"]
                pairs = pool[
                    rng.choice(len(pool), min(cfg["batch_size"] // 2, len(pool)), replace=False)
                ]
                loss = objective_loss(task, task["train"], pairs)
            else:
                pool = task["train"]
                batch = rng.choice(pool, min(len(pool), cfg["batch_size"]), replace=False)
                loss = objective_loss(task, batch)
            if not torch.isfinite(loss):
                raise ValueError(f"Non-finite training loss: {task['assay_id']}")
            optimizer.zero_grad(set_to_none=True)
            (loss * cfg[f"{objective}_weight"]).backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0, error_if_nonfinite=True)
            optimizer.step()
            train_losses[objective].append(float(loss.detach()))
        model.eval()
        prediction = np.full((len(val_ids), len(assays) * 2), np.nan, dtype=np.float32)
        with torch.inference_mode():
            for start in range(0, len(val_ids), cfg["batch_size"]):
                batch = val_ids[start : start + cfg["batch_size"]]
                prediction[start : start + len(batch)] = model(graph(batch)).cpu().numpy()
            validation_losses = {
                name: {source: [] for source in by_source[name]} for name in available
            }
            for task in tasks:
                positions, pairs = task["validation"], task["validation_pairs"]
                if not len(positions) or (pairs is not None and not len(pairs)):
                    continue
                value = float(objective_loss(task, positions, pairs, prediction))
                assert np.isfinite(value)
                validation_losses[task["objective"]][task["source"]].append(value)
        source_averages = {
            name: {source: float(np.mean(v)) for source, v in sources.items() if v}
            for name, sources in validation_losses.items()
        }
        averages = {
            name: float(np.mean(list(values.values())))
            for name, values in source_averages.items()
            if values
        }
        if not averages:
            raise ValueError("No held-out task losses for checkpoint selection")
        validation_loss = float(
            np.mean([v * cfg[f"{name}_weight"] for name, v in averages.items()])
        )
        row = dict(
            epoch=epoch,
            lr=lr,
            validation_loss=validation_loss,
            validation_by_objective=averages,
            validation_by_source=source_averages,
            train_by_objective={k: float(np.mean(v)) for k, v in train_losses.items() if v},
        )
        print(json.dumps(row), flush=True)
        with (out / "history.jsonl").open("a") as stream:
            stream.write(json.dumps(row) + "\n")
        checkpoint = dict(
            state_dict=model.state_dict(),
            optimizer=optimizer.state_dict(),
            epoch=epoch,
            config=cfg,
            assays=assays,
            tasks=task_meta,
            validation=row,
        )
        torch.save(checkpoint, out / "last.pt")
        if validation_loss < best:
            best, stale, best_epoch = validation_loss, 0, epoch
            best_state = copy.deepcopy(model.state_dict())
            torch.save(checkpoint, out / "best.pt")
        else:
            stale += 1
        if stale >= cfg["patience"]:
            break
    model.load_state_dict(best_state, strict=True)
    model.eval()
    values = []
    with torch.inference_mode():
        for start in range(0, challenge.height, cfg["batch_size"]):
            batch = BatchMolGraph(
                [
                    featurizer(Chem.MolFromSmiles(s))
                    for s in challenge.slice(start, cfg["batch_size"])["standardized_smiles"]
                ]
            )
            batch.to(device)
            values.append(model.embed(batch).cpu().numpy())
    values = np.concatenate(values)
    assert values.shape == (challenge.height, 600) and np.isfinite(values).all()
    pl.DataFrame(values, schema=[f"emb_{i:03d}" for i in range(600)]).insert_column(
        0, challenge["compound_id"]
    ).write_parquet(out / "external-embeddings.parquet")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                config=cfg,
                assays=len(assays),
                tasks=len(tasks),
                training_compounds=int((~is_val).sum()),
                validation_compounds=int(is_val.sum()),
                best_epoch=best_epoch,
                epochs=epoch + 1,
                best_validation_loss=best,
                challenge_rows=challenge.height,
                embedding_dimensions=600,
                challenge_labels_used=False,
                initialization="random",
                task_metadata=task_meta,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
