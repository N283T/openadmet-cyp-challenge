"""Train controlled Jeremy raw-embedding distillation students and audit transfer."""

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
import torch.nn.functional as functional
from chemprop import featurizers
from chemprop import nn as cnn
from chemprop.data import BatchMolGraph
from chemprop.nn.transforms import GraphTransform, UnscaleTransform
from rdkit import Chem

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def input_path(name: str) -> Path:
    path = inp / name
    if not path.is_file():
        raise ValueError(f"Missing declared input: {name}")
    return path


def read(name: str) -> pl.DataFrame:
    return pl.read_parquet(input_path(name))


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

    def forward(self, graph: BatchMolGraph) -> tuple[torch.Tensor, torch.Tensor]:
        embedding = self.embed(graph)
        return embedding, self.predictor(embedding)


def load_jeremy_checkpoint() -> dict:
    torch.serialization.add_safe_globals(
        [
            (dict, "lightning.fabric.utilities.data.AttributeDict"),
            GraphTransform,
            UnscaleTransform,
            cnn.BondMessagePassing,
            cnn.NormAggregation,
            cnn.MeanAggregation,
            cnn.RegressionFFN,
            torch.nn.Identity,
            torch.nn.ReLU,
            set,
        ]
    )
    checkpoint = torch.load(
        input_path("jeremy-checkpoint.pt"), map_location="cpu", weights_only=True
    )
    hyper = checkpoint["hyper_parameters"]
    assert hyper["message_passing"]["d_h"] == cfg["width"]
    assert hyper["message_passing"]["depth"] == cfg["depth"]
    assert hyper["agg"]["norm"] == cfg["aggregation_norm"]
    return checkpoint


def message_passing_state(checkpoint: dict) -> dict:
    return {
        key.removeprefix("message_passing."): value
        for key, value in checkpoint["state_dict"].items()
        if key.startswith("message_passing.")
    }


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
    result = []
    for smiles in frame["standardized_smiles"]:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError(f"Invalid standardized SMILES: {smiles}")
        result.append(featurizer(molecule))
    return result


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


def build_targets() -> tuple[np.ndarray, np.ndarray, np.ndarray, list[dict]]:
    observations = read("observations-compact-all.parquet").filter(
        pl.col("measurement_type") != "binary"
    )
    if task_sources := cfg.get("task_sources"):
        observations = observations.filter(pl.col("source").is_in(task_sources))
    observations = observations.with_columns(
        pl.col("value").cast(pl.Float64).alias("target")
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
        collapsed.with_columns(
            pl.col("row").is_in(np.flatnonzero(is_validation)).alias("is_validation")
        )
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
    selected = collapsed.filter(pl.col("assay_id").is_in(tasks))
    matrix = np.full((external.height, len(tasks)), np.nan, dtype=np.float32)
    task_index = {task: index for index, task in enumerate(tasks)}
    for row in selected.select("row", "assay_id", "target").iter_rows(named=True):
        matrix[row["row"], task_index[row["assay_id"]]] = row["target"]
    train_matrix = matrix[~is_validation]
    means = np.nanmean(train_matrix, axis=0).astype(np.float32)
    scales = np.nanstd(train_matrix, axis=0).astype(np.float32)
    if np.any(scales < 1e-8) or not np.isfinite(means).all() or not np.isfinite(scales).all():
        raise ValueError("Invalid task scaling")
    matrix = (matrix - means) / scales
    metadata = counts.with_columns(
        pl.Series("mean", means), pl.Series("scale", scales)
    ).to_dicts()
    return matrix, means, scales, metadata


targets, target_means, target_scales, task_metadata = build_targets()
pl.DataFrame(task_metadata).write_csv(out / "tasks.csv")
owned_checkpoint = None
if any(variant["initialization"] == "owned" for variant in cfg["variants"]):
    owned_checkpoint = torch.load(
        input_path("owned-checkpoint.pt"), map_location="cpu", weights_only=False
    )
    if owned_checkpoint["variant"]["name"] != "joint_all_no_binary":
        raise ValueError("Unexpected owned checkpoint variant")
    if [row["assay_id"] for row in owned_checkpoint["tasks"]] != [
        row["assay_id"] for row in task_metadata
    ]:
        raise ValueError("Owned checkpoint task order differs from current matrix")
jeremy_checkpoint = load_jeremy_checkpoint()
jeremy_state = message_passing_state(jeremy_checkpoint)


def embed_values(model: torch.nn.Module, graphs: list, *, use_embed: bool) -> np.ndarray:
    model.eval()
    chunks = []
    with torch.inference_mode():
        for start in range(0, len(graphs), cfg["batch_size"]):
            indices = np.arange(start, min(start + cfg["batch_size"], len(graphs)))
            batch = graph_batch(graphs, indices)
            values = model.embed(batch) if use_embed else model(batch)
            chunks.append(values.cpu().numpy())
    return np.concatenate(chunks).astype(np.float32)


class FrozenJeremy(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.message_passing = cnn.BondMessagePassing(
            d_h=cfg["width"], depth=cfg["depth"], dropout=0.0, bias=False
        )
        self.message_passing.load_state_dict(jeremy_state, strict=True)
        self.aggregation = cnn.NormAggregation(norm=cfg["aggregation_norm"])

    def forward(self, graph: BatchMolGraph) -> torch.Tensor:
        return self.aggregation(self.message_passing(graph), graph.batch)


print(json.dumps({"stage": "teacher_embeddings"}), flush=True)
teacher = FrozenJeremy().to(device)
teacher_external = embed_values(teacher, external_graphs, use_embed=False)
teacher_challenge = embed_values(teacher, challenge_graphs, use_embed=False)
del teacher
torch.cuda.empty_cache()
teacher_mean = teacher_external[~is_validation].mean(axis=0)
teacher_std = teacher_external[~is_validation].std(axis=0)
active_dimensions = teacher_std > 1e-5
if int(active_dimensions.sum()) < 500:
    raise ValueError("Unexpectedly few active Jeremy dimensions")
teacher_scale = np.where(active_dimensions, teacher_std, 1.0).astype(np.float32)
teacher_mean = teacher_mean.astype(np.float32)


def distillation_loss(
    student: torch.Tensor, teacher_values: torch.Tensor
) -> tuple[torch.Tensor, torch.Tensor]:
    active = torch.as_tensor(active_dimensions, device=device)
    mean = torch.as_tensor(teacher_mean, device=device)
    scale = torch.as_tensor(teacher_scale, device=device)
    student_z = ((student - mean) / scale)[:, active]
    teacher_z = ((teacher_values - mean) / scale)[:, active]
    coordinate = functional.mse_loss(student_z, teacher_z)
    student_unit = functional.normalize(student_z, dim=1)
    teacher_unit = functional.normalize(teacher_z, dim=1)
    geometry = functional.mse_loss(
        student_unit @ student_unit.T, teacher_unit @ teacher_unit.T
    )
    return coordinate, geometry


def learning_rate(step: int, total_steps: int, warmup_steps: int, maximum: float) -> float:
    if step < warmup_steps:
        fraction = step / max(1, warmup_steps)
        return cfg["init_lr"] + fraction * (maximum - cfg["init_lr"])
    fraction = (step - warmup_steps) / max(1, total_steps - warmup_steps)
    return maximum * (cfg["final_lr"] / maximum) ** fraction


def export_embeddings(model: JointModel, path: Path) -> tuple[pl.DataFrame, dict]:
    values = embed_values(model, challenge_graphs, use_embed=True)
    frame = pl.DataFrame(values, schema=[f"emb_{index:03d}" for index in range(values.shape[1])])
    frame.insert_column(0, challenge["compound_id"])
    frame.write_parquet(path)
    std = values.astype(np.float64).std(axis=0)
    centered = values.astype(np.float64) - values.astype(np.float64).mean(axis=0, keepdims=True)
    singular = np.linalg.svd(centered, compute_uv=False, full_matrices=False)
    proportions = singular**2 / np.sum(singular**2)
    entropy = -np.sum(proportions[proportions > 0] * np.log(proportions[proportions > 0]))
    active = active_dimensions
    normalized_error = np.mean(
        ((values[:, active] - teacher_challenge[:, active]) / teacher_scale[active]) ** 2
    )
    cosine = np.sum(values * teacher_challenge, axis=1) / np.maximum(
        np.linalg.norm(values, axis=1) * np.linalg.norm(teacher_challenge, axis=1), 1e-12
    )
    return frame, {
        "dead_dimensions": int(np.sum(std < 1e-8)),
        "effective_rank": float(np.exp(entropy)),
        "teacher_normalized_mse": float(normalized_error),
        "teacher_row_cosine": float(np.mean(cosine)),
    }


spec = importlib.util.spec_from_file_location("official", input_path("custom_scoring_functions.py"))
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
        endpoint_rows = []
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
            endpoint_rows.append(
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
        endpoint_frame = pl.concat(endpoint_rows)
        rows.extend(endpoint_rows)
        metrics.append(
            {
                "method": method,
                "endpoint": endpoint,
                "ST_RAE": float(
                    scorer.rae_soft_threshold_absolute_error(
                        endpoint_frame["y_true"].to_numpy(),
                        endpoint_frame["prediction"].to_numpy(),
                        y_true_lower=endpoint_frame["conf_low"].to_numpy(),
                        y_true_upper=endpoint_frame["conf_high"].to_numpy(),
                    )
                ),
                "n": endpoint_frame.height,
            }
        )
    return pl.concat(rows), metrics


def reference_metrics(alias: str, filename: str) -> list[dict]:
    frame = read(filename)
    result = []
    for endpoint in endpoints:
        rows = frame.filter(pl.col("endpoint") == endpoint)
        result.append(
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
    return result


train_indices = np.flatnonzero(~is_validation)
validation_indices = np.flatnonzero(is_validation)
results = []
all_metrics = []
for variant_index, variant in enumerate(cfg["variants"]):
    name = variant["name"]
    variant_dir = out / f"trial-{variant_index:02d}-{name}"
    variant_dir.mkdir()
    training_seed = cfg["seed"] + variant_index
    torch.manual_seed(training_seed)
    model = JointModel(targets.shape[1])
    if variant["initialization"] == "owned":
        assert owned_checkpoint is not None
        model.load_state_dict(owned_checkpoint["state_dict"], strict=True)
    elif variant["initialization"] == "jeremy":
        model.message_passing.load_state_dict(jeremy_state, strict=True)
    elif variant["initialization"] != "random":
        raise ValueError(f"Unknown initialization: {variant['initialization']}")
    model.to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=cfg["init_lr"])
    batches_per_epoch = math.ceil(len(train_indices) / cfg["batch_size"])
    total_steps = cfg["max_epochs"] * batches_per_epoch
    warmup_steps = cfg["warmup_epochs"] * batches_per_epoch
    rng = np.random.default_rng(training_seed)
    best_loss = float("inf")
    best_epoch = -1
    best_state = None
    stale = 0
    step = 0
    history = []
    for epoch in range(cfg["max_epochs"]):
        model.train()
        sums = dict(total=0.0, task=0.0, embedding=0.0, geometry=0.0)
        seen = 0
        shuffled = rng.permutation(train_indices)
        for start in range(0, len(train_indices), cfg["batch_size"]):
            batch_indices = shuffled[start : start + cfg["batch_size"]]
            if len(batch_indices) == 0:
                continue
            embedding, prediction = model(graph_batch(external_graphs, batch_indices))
            teacher_values = torch.as_tensor(teacher_external[batch_indices], device=device)
            embedding_loss, geometry_loss = distillation_loss(embedding, teacher_values)
            y = torch.as_tensor(targets[batch_indices], device=device)
            mask = torch.isfinite(y)
            task_loss = (
                functional.mse_loss(prediction[mask], y[mask])
                if mask.any()
                else torch.zeros((), device=device)
            )
            loss = (
                variant["task_weight"] * task_loss
                + variant["embedding_weight"] * embedding_loss
                + variant["geometry_weight"] * geometry_loss
            )
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0, error_if_nonfinite=True)
            optimizer.step()
            step += 1
            lr = learning_rate(step, total_steps, warmup_steps, variant["max_lr"])
            for group in optimizer.param_groups:
                group["lr"] = lr
            count = len(batch_indices)
            seen += count
            for key, value in [
                ("total", loss),
                ("task", task_loss),
                ("embedding", embedding_loss),
                ("geometry", geometry_loss),
            ]:
                sums[key] += float(value.detach()) * count
        model.eval()
        validation_sums = dict(total=0.0, task=0.0, embedding=0.0, geometry=0.0)
        validation_seen = 0
        with torch.inference_mode():
            for start in range(0, len(validation_indices), cfg["batch_size"]):
                indices = validation_indices[start : start + cfg["batch_size"]]
                embedding, prediction = model(graph_batch(external_graphs, indices))
                teacher_values = torch.as_tensor(teacher_external[indices], device=device)
                embedding_loss, geometry_loss = distillation_loss(embedding, teacher_values)
                y = torch.as_tensor(targets[indices], device=device)
                mask = torch.isfinite(y)
                task_loss = functional.mse_loss(prediction[mask], y[mask])
                loss = (
                    variant["task_weight"] * task_loss
                    + variant["embedding_weight"] * embedding_loss
                    + variant["geometry_weight"] * geometry_loss
                )
                count = len(indices)
                validation_seen += count
                for key, value in [
                    ("total", loss),
                    ("task", task_loss),
                    ("embedding", embedding_loss),
                    ("geometry", geometry_loss),
                ]:
                    validation_sums[key] += float(value) * count
        record = {
            "epoch": epoch,
            "learning_rate": optimizer.param_groups[0]["lr"],
            **{f"train_{key}": value / seen for key, value in sums.items()},
            **{
                f"validation_{key}": value / validation_seen
                for key, value in validation_sums.items()
            },
        }
        history.append(record)
        print(json.dumps({"variant": name, **record}), flush=True)
        if record["validation_total"] < best_loss - 1e-6:
            best_loss = record["validation_total"]
            best_epoch = epoch
            best_state = copy.deepcopy(model.state_dict())
            stale = 0
        else:
            stale += 1
            if stale >= cfg["patience"]:
                break
    if best_state is None:
        raise RuntimeError(f"No checkpoint selected for {name}")
    model.load_state_dict(best_state)
    torch.save(
        {
            "state_dict": best_state,
            "config": cfg,
            "variant": variant,
            "tasks": task_metadata,
            "target_means": target_means,
            "target_scales": target_scales,
            "best_epoch": best_epoch,
            "best_validation_loss": best_loss,
        },
        variant_dir / "checkpoint.pt",
    )
    (variant_dir / "history.json").write_text(json.dumps(history, indent=2) + "\n")
    embeddings, geometry = export_embeddings(model, variant_dir / "external-embeddings.parquet")
    oof, metrics = downstream_oof(embeddings, name)
    oof.write_parquet(variant_dir / "external-oof.parquet")
    result = {
        "variant": name,
        "initialization": variant["initialization"],
        "tasks": targets.shape[1],
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
    del model, optimizer, best_state
    torch.cuda.empty_cache()

all_metrics.extend(reference_metrics("owned_control", "owned-oof.parquet"))
all_metrics.extend(reference_metrics("jeremy_control", "jeremy-oof.parquet"))
pl.DataFrame(all_metrics).write_csv(out / "metrics.csv")
(out / "summary.json").write_text(
    json.dumps(
        {
            "config": cfg,
            "teacher_active_dimensions": int(active_dimensions.sum()),
            "results": results,
        },
        indent=2,
    )
    + "\n"
)
