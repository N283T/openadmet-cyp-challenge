"""Train project-owned joint sparse-matrix encoders and audit raw LGBM transfer."""

import copy
import importlib.util
import json
import math
import os
from pathlib import Path

import lightgbm as lgb
import numpy as np
import polars as pl
import torch
from chemprop import featurizers
from chemprop import nn as cnn
from chemprop.data import BatchMolGraph
from rdkit import Chem

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(name: str) -> pl.DataFrame:
    exact = inp / name
    paths = (
        [exact]
        if exact.exists()
        else [path for path in inp.iterdir() if path.name.partition("-")[2] == name]
    )
    if len(paths) != 1:
        raise ValueError(f"Expected one input for {name}: {paths}")
    return pl.read_parquet(paths[0])


class JointModel(torch.nn.Module):
    def __init__(self, n_tasks: int):
        super().__init__()
        self.message_passing = cnn.BondMessagePassing(
            d_h=cfg["width"], depth=cfg["depth"], dropout=0.0, bias=False
        )
        self.aggregation = cnn.NormAggregation(norm=cfg["aggregation_norm"])
        layers = []
        for _ in range(cfg["ffn_hidden_layers"]):
            layers.extend(
                [
                    torch.nn.Linear(cfg["width"], cfg["width"]),
                    torch.nn.ReLU(),
                    torch.nn.Dropout(cfg["dropout"]),
                ]
            )
        layers.append(torch.nn.Linear(cfg["width"], n_tasks))
        self.predictor = torch.nn.Sequential(*layers)

    def embed(self, graph: BatchMolGraph) -> torch.Tensor:
        return self.aggregation(self.message_passing(graph), graph.batch)

    def forward(self, graph: BatchMolGraph) -> torch.Tensor:
        return self.predictor(self.embed(graph))


torch.manual_seed(cfg["seed"])
np.random.seed(cfg["seed"])
torch.set_num_threads(cfg["num_threads"])
torch.backends.cuda.matmul.allow_tf32 = False
device = cfg["device"]

external = read("external-compounds.parquet").sort("compound_id")
challenge = read("compounds.parquet").sort("compound_id")
challenge_connectivity = read("challenge-connectivity.parquet")
if external["compound_id"].n_unique() != external.height:
    raise ValueError("External compound IDs are not unique")
if challenge["compound_id"].n_unique() != challenge.height:
    raise ValueError("Challenge compound IDs are not unique")
if set(external["connectivity"]) & set(challenge_connectivity["connectivity"]):
    raise ValueError("External corpus overlaps challenge connectivity")

featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()


def make_graphs(frame: pl.DataFrame) -> list:
    graphs = []
    for smiles in frame["standardized_smiles"]:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError(f"Invalid standardized SMILES: {smiles}")
        graphs.append(featurizer(molecule))
    return graphs


print(
    json.dumps(
        {"stage": "featurize", "external": external.height, "challenge": challenge.height}
    ),
    flush=True,
)
external_graphs = make_graphs(external)
challenge_graphs = make_graphs(challenge)


def graph_batch(graphs: list, indices: np.ndarray) -> BatchMolGraph:
    result = BatchMolGraph([graphs[int(index)] for index in indices])
    result.to(device)
    return result


split_rng = np.random.default_rng(cfg["validation_seed"])
connectivities = split_rng.permutation(sorted(external["connectivity"].unique()))
validation_connectivities = set(
    connectivities[: max(1, int(len(connectivities) * cfg["validation_fraction"]))]
)
is_validation = np.array(
    [value in validation_connectivities for value in external["connectivity"]], dtype=bool
)
external.select("compound_id", "connectivity").with_columns(
    pl.Series("role", np.where(is_validation, "validation", "training"))
).write_parquet(out / "membership.parquet")

compound_index = {compound_id: index for index, compound_id in enumerate(external["compound_id"])}


def build_targets(variant: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[dict]]:
    observation_file = variant["observation_file"]
    observations = read(observation_file)
    if sources := variant.get("sources"):
        observations = observations.filter(pl.col("source").is_in(sources))
    if excluded_sources := variant.get("exclude_sources"):
        observations = observations.filter(~pl.col("source").is_in(excluded_sources))
    if excluded_measurements := variant.get("exclude_measurement_types"):
        observations = observations.filter(
            ~pl.col("measurement_type").is_in(excluded_measurements)
        )
    if variant.get("cyp_only"):
        observations = observations.filter(pl.col("isoform").str.starts_with("CYP"))
    if not set(observations["compound_id"]).issubset(set(external["compound_id"])):
        raise ValueError(f"Unknown external compound in {variant['name']}")
    observations = observations.with_columns(
        pl.when(pl.col("binary_label").is_not_null())
        .then(pl.col("binary_label"))
        .otherwise(pl.col("value"))
        .cast(pl.Float64)
        .alias("target")
    ).filter(pl.col("target").is_finite())
    collapsed = observations.group_by("compound_id", "assay_id").agg(
        pl.col("target").median(),
        pl.col("source").first(),
        pl.col("isoform").first(),
        pl.col("measurement_type").first(),
    )
    collapsed = collapsed.with_columns(
        pl.col("compound_id").replace_strict(compound_index, return_dtype=pl.Int64).alias("row")
    )
    counts = (
        collapsed.with_columns(pl.col("row").is_in(np.flatnonzero(is_validation)).alias("is_validation"))
        .group_by("assay_id")
        .agg(
            (~pl.col("is_validation")).sum().alias("train_n"),
            pl.col("is_validation").sum().alias("validation_n"),
            pl.col("source").first(),
            pl.col("isoform").first(),
            pl.col("measurement_type").first(),
        )
        .filter(
            (pl.col("train_n") >= cfg["min_train_observations"])
            & (pl.col("validation_n") >= cfg["min_validation_observations"])
        )
        .sort("assay_id")
    )
    tasks = counts["assay_id"].to_list()
    if not tasks:
        raise ValueError(f"No eligible tasks for {observation_file}")
    selected = collapsed.filter(pl.col("assay_id").is_in(tasks))
    matrix = np.full((external.height, len(tasks)), np.nan, dtype=np.float32)
    task_index = {task: index for index, task in enumerate(tasks)}
    for row in selected.select("row", "assay_id", "target").iter_rows(named=True):
        matrix[row["row"], task_index[row["assay_id"]]] = row["target"]
    train_matrix = matrix[~is_validation]
    means = np.nanmean(train_matrix, axis=0).astype(np.float32)
    scales = np.nanstd(train_matrix, axis=0).astype(np.float32)
    if not np.isfinite(means).all() or not np.isfinite(scales).all() or np.any(scales < 1e-8):
        raise ValueError(f"Invalid task scaling for {observation_file}")
    matrix = (matrix - means) / scales
    metadata = counts.with_columns(
        pl.Series("mean", means), pl.Series("scale", scales)
    ).to_dicts()
    return matrix, means, scales, metadata


def learning_rate(step: int, total_steps: int, warmup_steps: int) -> float:
    if step < warmup_steps:
        fraction = step / max(1, warmup_steps)
        return cfg["init_lr"] + fraction * (cfg["max_lr"] - cfg["init_lr"])
    fraction = (step - warmup_steps) / max(1, total_steps - warmup_steps)
    return cfg["max_lr"] * (cfg["final_lr"] / cfg["max_lr"]) ** fraction


def export_embeddings(model: JointModel, path: Path) -> dict:
    model.eval()
    chunks = []
    with torch.inference_mode():
        for start in range(0, challenge.height, cfg["batch_size"]):
            indices = np.arange(start, min(start + cfg["batch_size"], challenge.height))
            chunks.append(model.embed(graph_batch(challenge_graphs, indices)).cpu().numpy())
    values = np.concatenate(chunks).astype(np.float32)
    embedding_frame = pl.DataFrame(
        values, schema=[f"emb_{index:03d}" for index in range(values.shape[1])]
    )
    embedding_frame.insert_column(0, challenge["compound_id"])
    embedding_frame.write_parquet(path)
    std = values.astype(np.float64).std(axis=0)
    centered = values.astype(np.float64) - values.astype(np.float64).mean(axis=0, keepdims=True)
    singular = np.linalg.svd(centered, compute_uv=False, full_matrices=False)
    proportions = singular**2 / np.sum(singular**2)
    entropy = -np.sum(proportions[proportions > 0] * np.log(proportions[proportions > 0]))
    return {
        "dead_dimensions": int(np.sum(std < 1e-8)),
        "median_feature_std": float(np.median(std)),
        "effective_rank": float(np.exp(entropy)),
        "participation_ratio": float(1 / np.sum(proportions**2)),
        "zero_fraction": float(np.mean(values == 0)),
        "median_active_dimensions": float(np.median(np.sum(values != 0, axis=1))),
    }


spec = importlib.util.spec_from_file_location("official", inp / "custom_scoring_functions.py")
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)
assignments = read("parent-assignments.parquet")
measurements = read("measurements.parquet").filter(
    pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition")
)
endpoints = sorted(measurements["endpoint"].unique())


def downstream_oof(embeddings: pl.DataFrame, method: str) -> tuple[pl.DataFrame, list[dict]]:
    columns = sorted(column for column in embeddings.columns if column.startswith("emb_"))
    rows = []
    metrics = []
    for endpoint in endpoints:
        table = (
            measurements.filter(pl.col("endpoint") == endpoint)
            .join(assignments, on="compound_id", validate="1:1")
            .join(embeddings, on="compound_id", validate="1:1")
            .sort("compound_id")
        )
        x = table.select(columns).to_numpy()
        y = table["value"].to_numpy()
        for fold in range(5):
            validation = (table["fold"].to_numpy() == fold) & ~table["is_anchor"].to_numpy()
            model = lgb.LGBMRegressor(
                objective="regression_l1",
                n_estimators=300,
                num_leaves=15,
                min_child_samples=20,
                learning_rate=0.05,
                reg_lambda=1.0,
                n_jobs=cfg["num_threads"],
                random_state=42 + fold,
                verbosity=-1,
                deterministic=True,
                force_col_wise=True,
            )
            model.fit(x[~validation], y[~validation])
            prediction = model.predict(x[validation]).astype(float)
            selected = table.filter(pl.Series(validation))
            rows.append(
                selected.select(
                    "compound_id",
                    "endpoint",
                    "fold",
                    "is_neighbor",
                    "conf_low",
                    "conf_high",
                    pl.col("value").alias("y_true"),
                ).with_columns(
                    pl.lit(method).alias("method"),
                    pl.lit("20260911T104703655337Z-a3923255").alias("split_id"),
                    pl.Series("prediction", prediction),
                )
            )
        endpoint_rows = pl.concat(rows).filter(pl.col("endpoint") == endpoint)
        metrics.append(
            {
                "method": method,
                "endpoint": endpoint,
                "ST_RAE": float(
                    scorer.rae_soft_threshold_absolute_error(
                        endpoint_rows["y_true"].to_numpy(),
                        endpoint_rows["prediction"].to_numpy(),
                        y_true_lower=endpoint_rows["conf_low"].to_numpy(),
                        y_true_upper=endpoint_rows["conf_high"].to_numpy(),
                    )
                ),
                "n": endpoint_rows.height,
            }
        )
    return pl.concat(rows), metrics


results = []
all_metrics = []
for variant_index, variant in enumerate(cfg["variants"]):
    variant_dir = out / f"trial-{variant_index:02d}-{variant['name']}"
    variant_dir.mkdir()
    targets, means, scales, task_metadata = build_targets(variant)
    task_weights = np.ones(targets.shape[1], dtype=np.float32)
    if focus_isoform := variant.get("focus_isoform"):
        focus_multiplier = float(variant["focus_multiplier"])
        if focus_multiplier <= 0:
            raise ValueError("focus_multiplier must be positive")
        focused = np.array(
            [row["isoform"] == focus_isoform for row in task_metadata], dtype=bool
        )
        if not focused.any():
            raise ValueError(f"No tasks found for focus isoform {focus_isoform}")
        task_weights[focused] = focus_multiplier
    for row, weight in zip(task_metadata, task_weights, strict=True):
        row["loss_weight"] = float(weight)
    pl.DataFrame(task_metadata).write_csv(variant_dir / "tasks.csv")
    labeled = np.flatnonzero(np.isfinite(targets).any(axis=1))
    train_indices = labeled[~is_validation[labeled]]
    validation_indices = labeled[is_validation[labeled]]
    training_seed = int(variant.get("training_seed", cfg["seed"] + variant_index))
    torch.manual_seed(training_seed)
    model = JointModel(targets.shape[1]).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=cfg["init_lr"])
    batches_per_epoch = math.ceil(len(train_indices) / cfg["batch_size"])
    total_steps = cfg["max_epochs"] * batches_per_epoch
    warmup_steps = cfg["warmup_epochs"] * batches_per_epoch
    rng = np.random.default_rng(training_seed)
    loss_weights = torch.as_tensor(task_weights, device=device).unsqueeze(0)
    best_loss = float("inf")
    best_epoch = -1
    stale = 0
    best_state = None
    history = []
    step = 0
    for epoch in range(cfg["max_epochs"]):
        model.train()
        shuffled = rng.permutation(train_indices)
        train_loss_sum = 0.0
        train_cells = 0
        for start in range(0, len(shuffled), cfg["batch_size"]):
            indices = shuffled[start : start + cfg["batch_size"]]
            y = torch.as_tensor(targets[indices], device=device)
            mask = torch.isfinite(y)
            prediction = model(graph_batch(external_graphs, indices))
            cell_weights = loss_weights.expand_as(y)[mask]
            loss = (
                torch.square(prediction[mask] - y[mask]) * cell_weights
            ).sum() / cell_weights.sum()
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0, error_if_nonfinite=True)
            optimizer.step()
            step += 1
            lr = learning_rate(step, total_steps, warmup_steps)
            for group in optimizer.param_groups:
                group["lr"] = lr
            cells = int(mask.sum())
            train_loss_sum += float(loss.detach()) * cells
            train_cells += cells
        model.eval()
        validation_loss_sum = 0.0
        validation_cells = 0
        with torch.inference_mode():
            for start in range(0, len(validation_indices), cfg["batch_size"]):
                indices = validation_indices[start : start + cfg["batch_size"]]
                y = torch.as_tensor(targets[indices], device=device)
                mask = torch.isfinite(y)
                prediction = model(graph_batch(external_graphs, indices))
                cell_weights = loss_weights.expand_as(y)[mask]
                loss_sum = (
                    torch.square(prediction[mask] - y[mask]) * cell_weights
                ).sum()
                validation_loss_sum += float(loss_sum)
                validation_cells += float(cell_weights.sum())
        validation_loss = validation_loss_sum / validation_cells
        record = {
            "epoch": epoch,
            "learning_rate": optimizer.param_groups[0]["lr"],
            "train_loss": train_loss_sum / train_cells,
            "validation_loss": validation_loss,
        }
        history.append(record)
        print(json.dumps({"variant": variant["name"], **record}), flush=True)
        if validation_loss < best_loss - 1e-6:
            best_loss = validation_loss
            best_epoch = epoch
            stale = 0
            best_state = copy.deepcopy(model.state_dict())
        else:
            stale += 1
            if stale >= cfg["patience"]:
                break
    if best_state is None:
        raise RuntimeError(f"No checkpoint selected for {variant['name']}")
    model.load_state_dict(best_state)
    torch.save(
        {
            "state_dict": best_state,
            "config": cfg,
            "variant": variant,
            "tasks": task_metadata,
            "target_means": means,
            "target_scales": scales,
            "best_epoch": best_epoch,
            "best_validation_loss": best_loss,
        },
        variant_dir / "checkpoint.pt",
    )
    (variant_dir / "history.json").write_text(json.dumps(history, indent=2) + "\n")
    geometry = export_embeddings(model, variant_dir / "external-embeddings.parquet")
    embeddings = pl.read_parquet(variant_dir / "external-embeddings.parquet")
    oof, metrics = downstream_oof(embeddings, variant["name"])
    oof.write_parquet(variant_dir / "external-oof.parquet")
    result = {
        "variant": variant["name"],
        "observation_file": variant["observation_file"],
        "tasks": targets.shape[1],
        "training_molecules": len(train_indices),
        "validation_molecules": len(validation_indices),
        "best_epoch": best_epoch,
        "best_validation_loss": best_loss,
        "geometry": geometry,
        "metrics": metrics,
        "macro_ST_RAE": float(np.mean([row["ST_RAE"] for row in metrics])),
    }
    results.append(result)
    all_metrics.extend(metrics)
    (out / "summary.json").write_text(
        json.dumps({"config": cfg, "results": results}, indent=2) + "\n"
    )
    print(json.dumps(result), flush=True)
    del model, optimizer, targets, best_state
    torch.cuda.empty_cache()


def reference_metrics(alias: str, name: str) -> list[dict]:
    frame = read(name)
    values = []
    for endpoint in endpoints:
        rows = frame.filter(pl.col("endpoint") == endpoint)
        values.append(
            {
                "method": alias,
                "endpoint": endpoint,
                "ST_RAE": float(
                    scorer.rae_soft_threshold_absolute_error(
                        rows["y_true"].to_numpy(),
                        rows["prediction"].to_numpy(),
                        y_true_lower=rows["conf_low"].to_numpy(),
                        y_true_upper=rows["conf_high"].to_numpy(),
                    )
                ),
                "n": rows.height,
            }
        )
    return values


all_metrics.extend(reference_metrics("previous_owned", "previous-owned.parquet"))
all_metrics.extend(reference_metrics("jeremy", "jeremy-oof.parquet"))
pl.DataFrame(all_metrics).write_csv(out / "metrics.csv")
(out / "summary.json").write_text(json.dumps({"config": cfg, "results": results}, indent=2) + "\n")
