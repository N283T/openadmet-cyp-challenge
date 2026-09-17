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
from external_assay_scaling import GROUP_COLUMNS, shrunk_assay_std, training_variance_pools
from external_loss_utils import interval_mse, ordered_pairs, ranking_loss
from rdkit import Chem


def read_input(root, name):
    exact = root / name
    paths = [exact] if exact.exists() else [
        p for p in root.iterdir() if p.name.partition("-")[2] == name
    ]
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
        scaling_policy="assay",
        scaling_pseudocount=20.0,
        scaling_minimum_std=0.3,
    )
    cfg.update(json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text()))
    if cfg["scaling_policy"] not in {"assay", "pooled_shrinkage"}:
        raise ValueError("Unknown scaling_policy")
    modes = {
        "continuous": ["continuous"],
        "continuous_ranking": ["continuous", "ranking"],
        "ranking": ["ranking"],
        "classification": ["classification"],
        "mixed": ["continuous", "ranking", "classification"],
    }
    objectives = modes[cfg["mode"]]
    (out / "effective-config.json").write_text(json.dumps(cfg, indent=2))
    compounds = read_input(inp, "external-compounds.parquet").sort("compound_id")
    observations = read_input(inp, "observations.parquet")
    task_filter = cfg.get("task_filter", {})
    if prefix := task_filter.get("isoform_prefix"):
        observations = observations.filter(pl.col("isoform").str.starts_with(prefix))
    if target_isoform := task_filter.get("target_cyp_isoform"):
        if not target_isoform.startswith("CYP"):
            raise ValueError(f"target_cyp_isoform must start with CYP: {target_isoform}")
        observations = observations.filter(
            ~pl.col("isoform").str.starts_with("CYP")
            | (pl.col("isoform") == target_isoform)
        )
    if sources := task_filter.get("sources"):
        observations = observations.filter(pl.col("source").is_in(sources))
    if measurement_types := task_filter.get("measurement_types"):
        observations = observations.filter(
            pl.col("measurement_type").is_in(measurement_types)
        )
    if not observations.height:
        raise ValueError(f"Task filter removed every observation: {task_filter}")
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
            *[
                pl.col(name).first()
                for name in ["isoform", "measurement_type"]
                if name in observations.columns
            ],
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
    validation_seed = cfg.get("validation_seed", cfg["seed"] + 739)
    val_rng = np.random.default_rng(validation_seed)
    split_rng = np.random.default_rng(validation_seed)
    keys = split_rng.permutation(sorted(compounds["connectivity"].unique()))
    val_keys = set(keys[: max(1, int(len(keys) * cfg["validation_fraction"]))])
    is_val = np.array([key in val_keys for key in compounds["connectivity"]])
    compounds.select("compound_id", "connectivity").with_columns(
        pl.Series("role", np.where(is_val, "validation", "training"))
    ).write_parquet(out / "membership.parquet")
    scaling_pools = {}
    if cfg["scaling_policy"] == "pooled_shrinkage":
        for column in GROUP_COLUMNS:
            if column not in observations.columns:
                raise ValueError(f"pooled_shrinkage requires {column}")
        if (
            observations.group_by("assay_id")
            .agg(pl.struct(GROUP_COLUMNS).n_unique().alias("metadata_variants"))
            .filter(pl.col("metadata_variants") != 1)
            .height
        ):
            raise ValueError("Assay metadata must map to one source/CYP/measurement group")
        pool_table = training_variance_pools(
            observations,
            compounds["compound_id"].filter(pl.Series(~is_val)).to_list(),
            cfg["scaling_minimum_std"],
        )
        pool_table.write_parquet(out / "training-variance-pools.parquet")
        scaling_pools = {
            tuple(row[column] for column in GROUP_COLUMNS): row
            for row in pool_table.iter_rows(named=True)
        }
    tasks, task_meta = [], []
    assays = sorted(observations["assay_id"].unique())
    assay_metadata = {
        row["assay_id"]: row
        for row in observations.group_by("assay_id").agg(
            pl.col("isoform").first(),
            pl.col("measurement_type").first(),
            pl.col("isoform").n_unique().alias("isoform_variants"),
            pl.col("measurement_type").n_unique().alias("measurement_type_variants"),
        ).iter_rows(named=True)
    }
    if any(
        row["isoform_variants"] != 1 or row["measurement_type_variants"] != 1
        for row in assay_metadata.values()
    ):
        raise ValueError("Each assay must have one isoform and measurement type")
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
        scaling_meta = dict(scaling_policy=cfg["scaling_policy"])
        if cfg["scaling_policy"] == "pooled_shrinkage":
            group = tuple(rows[column][0] for column in GROUP_COLUMNS)
            pool = scaling_pools.get(group)
            # Ranking/binary-only groups can have no exact training observations.
            pool_variance = pool["pool_variance"] if pool else 1.0
            std = shrunk_assay_std(
                lower[train_exact],
                pool_variance,
                cfg["scaling_pseudocount"],
                cfg["scaling_minimum_std"],
            )
            scaling_meta.update(
                scaling_group=dict(zip(GROUP_COLUMNS, group, strict=True)),
                scaling_pool_variance=pool_variance,
                scaling_pool_exact_n=pool["pool_exact_n"] if pool else 0,
                scaling_pool_fallback=pool is None,
                scaling_pseudocount=cfg["scaling_pseudocount"],
                scaling_minimum_std=cfg["scaling_minimum_std"],
                scaling_train_exact_n=int(train_exact.sum()),
            )
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
                    **scaling_meta,
                )
            )
    if not tasks:
        raise ValueError("No eligible training tasks")
    # Allocate outputs only for eligible tasks. Optionally let semantically matched
    # assays share a standardized output head while retaining assay-specific losses.
    old_assays = assays
    assays = sorted({task["assay_id"] for task in tasks})
    head_grouping = cfg.get("head_grouping", "assay")
    if head_grouping == "assay":
        head_ids = {assay: assay for assay in assays}
    elif head_grouping == "isoform_measurement":
        head_ids = {
            assay: (
                f"{assay_metadata[assay]['isoform']}::"
                f"{assay_metadata[assay]['measurement_type']}"
            )
            for assay in assays
        }
    else:
        raise ValueError(f"Unknown head_grouping: {head_grouping}")
    heads = sorted(set(head_ids.values()))
    head_lookup = {head: i for i, head in enumerate(heads)}
    for task, meta in zip(tasks, task_meta, strict=True):
        head_id = head_ids[task["assay_id"]]
        task["head_id"] = head_id
        meta["head_id"] = head_id
        task["column"] = head_lookup[head_id] * 2 + int(
            task["objective"] == "classification"
        )
        meta["column"] = task["column"]
    print(
        json.dumps(
            dict(
                raw_assays=len(old_assays),
                eligible_assays=len(assays),
                output_heads=len(heads),
                tasks=len(tasks),
            )
        ),
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

    model = AssayModel(len(heads), cfg["dropout"]).to(device)
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
    sampling_weights = np.array(
        [cfg.get("objective_sampling_weights", {}).get(name, 1.0) for name in available],
        dtype=float,
    )
    if not np.isfinite(sampling_weights).all() or np.any(sampling_weights <= 0):
        raise ValueError(f"Invalid objective sampling weights: {sampling_weights}")
    sampling_weights /= sampling_weights.sum()
    by_source = {
        name: {
            source: [t for t in ts if t["source"] == source]
            for source in sorted({t["source"] for t in ts})
        }
        for name, ts in available.items()
    }
    task_sampling = cfg.get("task_sampling", "source_uniform")
    if task_sampling not in {
        "source_uniform",
        "task_uniform",
        "sqrt_observations",
        "observation_power",
    }:
        raise ValueError(f"Unknown task_sampling: {task_sampling}")
    task_probabilities = {}
    for name, objective_tasks in available.items():
        if task_sampling == "task_uniform":
            weights = np.ones(len(objective_tasks), dtype=float)
        elif task_sampling == "sqrt_observations":
            weights = np.sqrt([len(task["train"]) for task in objective_tasks])
        elif task_sampling == "observation_power":
            power = float(cfg["task_sampling_power"])
            if not 0 <= power <= 1:
                raise ValueError(f"task_sampling_power must be in [0, 1]; got {power}")
            weights = np.power([len(task["train"]) for task in objective_tasks], power)
        else:
            continue
        task_probabilities[name] = weights / weights.sum()
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
            objective = rng.choice(list(available), p=sampling_weights)
            if task_sampling == "source_uniform":
                source = rng.choice(list(by_source[objective]))
                pool_tasks = by_source[objective][source]
                task = pool_tasks[rng.integers(len(pool_tasks))]
            else:
                task_index = rng.choice(
                    len(available[objective]), p=task_probabilities[objective]
                )
                task = available[objective][task_index]
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
        prediction = np.full((len(val_ids), len(heads) * 2), np.nan, dtype=np.float32)
        with torch.inference_mode():
            for start in range(0, len(val_ids), cfg["batch_size"]):
                batch = val_ids[start : start + cfg["batch_size"]]
                prediction[start : start + len(batch)] = model(graph(batch)).cpu().numpy()
            validation_losses = {
                name: {source: [] for source in by_source[name]} for name in available
            }
            validation_task_rows = []
            for task in tasks:
                positions, pairs = task["validation"], task["validation_pairs"]
                if not len(positions) or (pairs is not None and not len(pairs)):
                    continue
                value = float(objective_loss(task, positions, pairs, prediction))
                assert np.isfinite(value)
                validation_losses[task["objective"]][task["source"]].append(value)
                validation_task_rows.append(
                    {
                        "assay_id": task["assay_id"],
                        "source": task["source"],
                        "isoform": assay_metadata[task["assay_id"]]["isoform"],
                        "measurement_type": assay_metadata[task["assay_id"]][
                            "measurement_type"
                        ],
                        "objective": task["objective"],
                        "validation_n": len(positions),
                        "loss": value,
                    }
                )
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
        checkpoint_rows = validation_task_rows
        checkpoint_filter = cfg.get("checkpoint_task_filter", {})
        if names := checkpoint_filter.get("objectives"):
            checkpoint_rows = [row for row in checkpoint_rows if row["objective"] in names]
        if prefix := checkpoint_filter.get("isoform_prefix"):
            checkpoint_rows = [
                row for row in checkpoint_rows if row["isoform"].startswith(prefix)
            ]
        if sources := checkpoint_filter.get("sources"):
            checkpoint_rows = [row for row in checkpoint_rows if row["source"] in sources]
        if not checkpoint_rows:
            raise ValueError(f"No checkpoint tasks matched {checkpoint_filter}")
        checkpoint_aggregation = cfg.get("checkpoint_aggregation", "source_objective")
        if checkpoint_aggregation == "source_objective" and not checkpoint_filter:
            validation_loss = float(
                np.mean([v * cfg[f"{name}_weight"] for name, v in averages.items()])
            )
        elif checkpoint_aggregation == "task_macro":
            validation_loss = float(np.mean([row["loss"] for row in checkpoint_rows]))
        else:
            raise ValueError(
                "Filtered checkpoint tasks require checkpoint_aggregation=task_macro; "
                f"got {checkpoint_aggregation}"
            )
        validation_observation_weighted = float(
            np.average(
                [row["loss"] for row in checkpoint_rows],
                weights=[row["validation_n"] for row in checkpoint_rows],
            )
        )
        row = dict(
            epoch=epoch,
            lr=lr,
            validation_loss=validation_loss,
            validation_observation_weighted=validation_observation_weighted,
            checkpoint_task_count=len(checkpoint_rows),
            checkpoint_task_filter=checkpoint_filter,
            checkpoint_aggregation=checkpoint_aggregation,
            validation_by_objective=averages,
            validation_by_source=source_averages,
            validation_by_task=validation_task_rows,
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
            output_heads=heads,
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
    embedding_dimensions = values.shape[1]
    assert values.shape[0] == challenge.height and np.isfinite(values).all()
    pl.DataFrame(
        values, schema=[f"emb_{i:04d}" for i in range(embedding_dimensions)]
    ).insert_column(0, challenge["compound_id"]).write_parquet(out / "external-embeddings.parquet")
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                config=cfg,
                assays=len(assays),
                output_heads=len(heads),
                head_grouping=head_grouping,
                tasks=len(tasks),
                training_compounds=int((~is_val).sum()),
                validation_compounds=int(is_val.sum()),
                best_epoch=best_epoch,
                epochs=epoch + 1,
                best_validation_loss=best,
                challenge_rows=challenge.height,
                embedding_dimensions=embedding_dimensions,
                validation_seed=validation_seed,
                challenge_labels_used=False,
                initialization="random",
                task_metadata=task_meta,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
