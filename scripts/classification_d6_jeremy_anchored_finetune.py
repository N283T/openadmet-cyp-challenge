"""Fold-safe Jeremy-anchored D6 fine-tuning with CheMeleon readouts."""

import copy
import json
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
from sklearn.metrics import average_precision_score, matthews_corrcoef, roc_auc_score
from sklearn.model_selection import StratifiedShuffleSplit

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def read(name: str) -> pl.DataFrame:
    path = inp / name
    if not path.is_file():
        raise FileNotFoundError(path)
    return pl.read_parquet(path)


class EncoderClassifier(torch.nn.Module):
    def __init__(self, state: dict, width: int, depth: int, norm: float):
        super().__init__()
        self.message_passing = cnn.BondMessagePassing(
            d_h=width, depth=depth, dropout=0.0, bias=False
        )
        self.message_passing.load_state_dict(state, strict=True)
        self.aggregation = cnn.NormAggregation(norm=norm)
        self.head = torch.nn.Linear(width, 1)

    def embed(self, graph: BatchMolGraph) -> torch.Tensor:
        return self.aggregation(self.message_passing(graph), graph.batch)

    def forward(self, graph: BatchMolGraph) -> tuple[torch.Tensor, torch.Tensor]:
        embedding = self.embed(graph)
        return embedding, self.head(embedding).squeeze(1)


def load_jeremy() -> tuple[dict, dict]:
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
        inp / "jeremy-checkpoint.pt", map_location="cpu", weights_only=True
    )
    state = {
        key.removeprefix("message_passing."): value
        for key, value in checkpoint["state_dict"].items()
        if key.startswith("message_passing.")
    }
    return checkpoint["hyper_parameters"], state


def fixed_calls(scores: np.ndarray, folds: np.ndarray) -> np.ndarray:
    calls = np.zeros(len(scores), dtype=bool)
    for fold in np.unique(folds):
        positions = np.flatnonzero(folds == fold)
        count = max(1, int(round(cfg["call_rate"] * len(positions))))
        order = np.argsort(-scores[positions], kind="stable")
        calls[positions[order[:count]]] = True
    return calls


def expected_mcc(y: np.ndarray, calls: np.ndarray) -> float:
    positive = y == 1
    tpr = float(calls[positive].mean())
    fpr = float(calls[~positive].mean())
    prevalence = cfg["target_prior"]
    tp, fn = prevalence * tpr, prevalence * (1 - tpr)
    fp, tn = (1 - prevalence) * fpr, (1 - prevalence) * (1 - fpr)
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float((tp * tn - fp * fn) / denominator) if denominator else 0.0


torch.manual_seed(cfg["seed"])
np.random.seed(cfg["seed"])
torch.set_num_threads(cfg["num_threads"])
torch.backends.cuda.matmul.allow_tf32 = False
torch.backends.cudnn.allow_tf32 = False
device = torch.device(cfg["device"])

hyper, jeremy_state = load_jeremy()
width = int(hyper["message_passing"]["d_h"])
depth = int(hyper["message_passing"]["depth"])
norm = float(hyper["agg"]["norm"])
compounds = read("compounds.parquet").sort("compound_id")
canonical_ids = compounds["compound_id"].to_list()
id_to_index = {compound_id: index for index, compound_id in enumerate(canonical_ids)}
chemeleon = read("chemeleon-raw.parquet").sort("compound_id")
if chemeleon["compound_id"].to_list() != canonical_ids:
    raise ValueError("CheMeleon compound order mismatch")
chemeleon_columns = sorted(
    column for column in chemeleon.columns if column.startswith("emb_")
)
chemeleon_values = chemeleon.select(chemeleon_columns).to_numpy().astype(np.float32)

labels = (
    read("measurements.parquet")
    .filter(pl.col("endpoint") == cfg["endpoint"])
    .select("compound_id", pl.col("value").cast(pl.Int8).alias("y_true"))
    .join(
        read("assignments.parquet")
        .filter(pl.col("strategy") == cfg["strategy"])
        .select("compound_id", "fold", "group"),
        on="compound_id",
        validate="1:1",
    )
    .sort("compound_id")
)
label_indices = np.asarray([id_to_index[value] for value in labels["compound_id"]])
y = labels["y_true"].to_numpy().astype(np.float32)
folds = labels["fold"].to_numpy().astype(int)
test = read("test-mapping.parquet").select("compound_id").unique().sort("compound_id")
test_indices = np.asarray([id_to_index[value] for value in test["compound_id"]])

featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
graphs = []
for compound_id, smiles in compounds.select(
    "compound_id", "standardized_smiles"
).iter_rows():
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Invalid standardized SMILES: {compound_id}")
    graphs.append(featurizer(molecule))


def graph_batch(indices: np.ndarray) -> BatchMolGraph:
    graph = BatchMolGraph([graphs[int(index)] for index in indices])
    graph.to(device)
    return graph


teacher_model = EncoderClassifier(jeremy_state, width, depth, norm).to(device)
teacher_model.head = torch.nn.Identity()
teacher_model.eval()


def teacher_embeddings(model: EncoderClassifier, indices: np.ndarray) -> np.ndarray:
    chunks = []
    with torch.inference_mode():
        for start in range(0, len(indices), cfg["batch_size"]):
            batch = indices[start : start + cfg["batch_size"]]
            chunks.append(model.embed(graph_batch(batch)).cpu().numpy())
    return np.concatenate(chunks).astype(np.float32)


all_indices = np.arange(len(canonical_ids))
teacher_values = teacher_embeddings(teacher_model, all_indices)
teacher_mean = teacher_values.mean(axis=0).astype(np.float32)
teacher_scale = teacher_values.std(axis=0).astype(np.float32)
active = teacher_scale > 1e-5
teacher_scale = np.where(active, teacher_scale, 1.0).astype(np.float32)
del teacher_model
torch.cuda.empty_cache()


def anchor_losses(student: torch.Tensor, teacher: torch.Tensor):
    active_tensor = torch.as_tensor(active, device=device)
    mean = torch.as_tensor(teacher_mean, device=device)
    scale = torch.as_tensor(teacher_scale, device=device)
    student_z = ((student - mean) / scale)[:, active_tensor]
    teacher_z = ((teacher - mean) / scale)[:, active_tensor]
    coordinate = functional.mse_loss(student_z, teacher_z)
    student_unit = functional.normalize(student_z, dim=1)
    teacher_unit = functional.normalize(teacher_z, dim=1)
    geometry = functional.mse_loss(
        student_unit @ student_unit.T, teacher_unit @ teacher_unit.T
    )
    return coordinate, geometry


def embed_all(model: EncoderClassifier) -> np.ndarray:
    chunks = []
    model.eval()
    with torch.inference_mode():
        for start in range(0, len(all_indices), cfg["batch_size"]):
            indices = all_indices[start : start + cfg["batch_size"]]
            chunks.append(model.embed(graph_batch(indices)).cpu().numpy())
    return np.concatenate(chunks).astype(np.float32)


oof_scores = {
    f"{variant['name']}::{backend}": np.empty(len(labels), dtype=float)
    for variant in cfg["variants"]
    for backend in ("direct", "lgbm_balanced")
}
test_fold_scores = {method: [] for method in oof_scores}
history_rows = []
geometry_rows = []

for fold in range(5):
    outer_train = np.flatnonzero(folds != fold)
    outer_validation = np.flatnonzero(folds == fold)
    splitter = StratifiedShuffleSplit(
        n_splits=1,
        test_size=cfg["inner_validation_fraction"],
        random_state=cfg["seed"] + fold,
    )
    train_relative, inner_relative = next(
        splitter.split(np.zeros(len(outer_train)), y[outer_train])
    )
    train_rows = outer_train[train_relative]
    inner_rows = outer_train[inner_relative]
    for variant_index, variant in enumerate(cfg["variants"]):
        torch.manual_seed(cfg["seed"] + 100 * fold + variant_index)
        model = EncoderClassifier(jeremy_state, width, depth, norm).to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=cfg["learning_rate"])
        positive = float(y[train_rows].sum())
        pos_weight = torch.tensor(
            [(len(train_rows) - positive) / positive], device=device
        )
        best_loss = float("inf")
        best_state = None
        stale = 0
        rng = np.random.default_rng(cfg["seed"] + 100 * fold + variant_index)
        for epoch in range(cfg["max_epochs"]):
            model.train()
            train_losses = []
            shuffled_rows = rng.permutation(train_rows)
            for start in range(0, len(train_rows), cfg["batch_size"]):
                batch_rows = shuffled_rows[start : start + cfg["batch_size"]]
                graph_indices = label_indices[batch_rows]
                embedding, logits = model(graph_batch(graph_indices))
                target = torch.as_tensor(y[batch_rows], device=device)
                task_loss = functional.binary_cross_entropy_with_logits(
                    logits, target, pos_weight=pos_weight
                )
                teacher = torch.as_tensor(teacher_values[graph_indices], device=device)
                coordinate, geometry = anchor_losses(embedding, teacher)
                loss = (
                    task_loss
                    + variant["coordinate_weight"] * coordinate
                    + variant["geometry_weight"] * geometry
                )
                optimizer.zero_grad(set_to_none=True)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
                optimizer.step()
                train_losses.append(float(loss.detach()))
            model.eval()
            validation_losses = []
            with torch.inference_mode():
                for start in range(0, len(inner_rows), cfg["batch_size"]):
                    batch_rows = inner_rows[start : start + cfg["batch_size"]]
                    graph_indices = label_indices[batch_rows]
                    embedding, logits = model(graph_batch(graph_indices))
                    target = torch.as_tensor(y[batch_rows], device=device)
                    task_loss = functional.binary_cross_entropy_with_logits(
                        logits, target, pos_weight=pos_weight
                    )
                    teacher = torch.as_tensor(
                        teacher_values[graph_indices], device=device
                    )
                    coordinate, geometry = anchor_losses(embedding, teacher)
                    validation_losses.append(
                        float(
                            task_loss
                            + variant["coordinate_weight"] * coordinate
                            + variant["geometry_weight"] * geometry
                        )
                    )
            validation_loss = float(np.mean(validation_losses))
            history_rows.append(
                {
                    "variant": variant["name"],
                    "fold": fold,
                    "epoch": epoch,
                    "train_loss": float(np.mean(train_losses)),
                    "validation_loss": validation_loss,
                }
            )
            print(
                json.dumps(
                    {
                        "variant": variant["name"],
                        "fold": fold,
                        "epoch": epoch,
                        "validation_loss": validation_loss,
                    }
                ),
                flush=True,
            )
            if validation_loss < best_loss - 1e-5:
                best_loss = validation_loss
                best_state = copy.deepcopy(model.state_dict())
                stale = 0
            else:
                stale += 1
                if stale >= cfg["patience"]:
                    break
        if best_state is None:
            raise RuntimeError(f"No checkpoint selected for {variant['name']} fold {fold}")
        model.load_state_dict(best_state)
        adapted = embed_all(model)
        normalized_mse = float(
            np.mean(
                (
                    (adapted[:, active] - teacher_values[:, active])
                    / teacher_scale[active]
                )
                ** 2
            )
        )
        cosine = float(
            np.mean(
                np.sum(adapted * teacher_values, axis=1)
                / np.maximum(
                    np.linalg.norm(adapted, axis=1)
                    * np.linalg.norm(teacher_values, axis=1),
                    1e-12,
                )
            )
        )
        geometry_rows.append(
            {
                "variant": variant["name"],
                "fold": fold,
                "best_validation_loss": best_loss,
                "teacher_normalized_mse": normalized_mse,
                "teacher_row_cosine": cosine,
            }
        )
        model.eval()
        with torch.inference_mode():
            _, validation_logits = model(graph_batch(label_indices[outer_validation]))
            _, test_logits = model(graph_batch(test_indices))
        direct_method = f"{variant['name']}::direct"
        oof_scores[direct_method][outer_validation] = (
            validation_logits.sigmoid().cpu().numpy()
        )
        test_fold_scores[direct_method].append(test_logits.sigmoid().cpu().numpy())

        features = np.column_stack([adapted, chemeleon_values])
        params = cfg["lightgbm"] | {
            "objective": "binary",
            "class_weight": "balanced",
            "n_jobs": cfg["num_threads"],
            "random_state": cfg["seed"] + fold,
            "verbosity": -1,
            "deterministic": True,
            "force_col_wise": True,
        }
        readout = lgb.LGBMClassifier(**params)
        readout.fit(features[label_indices[outer_train]], y[outer_train])
        readout_method = f"{variant['name']}::lgbm_balanced"
        oof_scores[readout_method][outer_validation] = readout.predict_proba(
            features[label_indices[outer_validation]]
        )[:, 1]
        test_fold_scores[readout_method].append(
            readout.predict_proba(features[test_indices])[:, 1]
        )
        checkpoint_dir = out / "checkpoints"
        checkpoint_dir.mkdir(exist_ok=True)
        torch.save(
            {
                "state_dict": best_state,
                "variant": variant,
                "fold": fold,
                "best_validation_loss": best_loss,
            },
            checkpoint_dir / f"{variant['name']}-fold{fold}.pt",
        )
        del model, optimizer, best_state, adapted, features, readout
        torch.cuda.empty_cache()

oof_rows = []
test_rows = []
metric_rows = []
for method, scores in oof_scores.items():
    calls = fixed_calls(scores, folds)
    metric_rows.append(
        {
            "method": method,
            "AUROC": float(roc_auc_score(y, scores)),
            "AP": float(average_precision_score(y, scores)),
            "MCC_observed": float(matthews_corrcoef(y, calls)),
            "MCC_target_prior": expected_mcc(y, calls),
            "positive_calls": int(calls.sum()),
        }
    )
    for row, score in enumerate(scores):
        oof_rows.append(
            {
                "compound_id": labels["compound_id"][row],
                "endpoint": cfg["endpoint"],
                "split_id": cfg["strategy"],
                "fold": int(folds[row]),
                "group": labels["group"][row],
                "y_true": int(y[row]),
                "method": method,
                "prediction": float(score),
                "call": bool(calls[row]),
            }
        )
    test_scores = np.mean(test_fold_scores[method], axis=0)
    test_calls = np.zeros(len(test), dtype=bool)
    count = int(round(cfg["call_rate"] * len(test)))
    test_calls[np.argsort(-test_scores, kind="stable")[:count]] = True
    for row, score in enumerate(test_scores):
        test_rows.append(
            {
                "compound_id": test["compound_id"][row],
                "endpoint": cfg["endpoint"],
                "split_id": cfg["strategy"],
                "fold": -1,
                "method": method,
                "prediction": float(score),
                "call": bool(test_calls[row]),
            }
        )

metrics = pl.DataFrame(metric_rows).sort("MCC_target_prior", descending=True)
metrics.write_csv(out / "metrics.csv")
pl.DataFrame(oof_rows).write_parquet(out / "oof.parquet")
pl.DataFrame(test_rows).write_parquet(out / "test.parquet")
pl.DataFrame(history_rows).write_parquet(out / "history.parquet")
pl.DataFrame(geometry_rows).write_csv(out / "geometry.csv")
(out / "summary.json").write_text(
    json.dumps({"config": cfg, "metrics": metrics.to_dicts()}, indent=2) + "\n"
)
print(metrics, flush=True)
