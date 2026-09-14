"""Train PXR Chemprop without any CYP challenge connectivity."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
from chemprop import data as chemprop_data
from chemprop import models, nn
from chemprop.nn.metrics import MSE
from lightning import pytorch as lightning
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
from scipy.stats import spearmanr
from sklearn.model_selection import GroupShuffleSplit

from cyp_challenge.preprocessing import standardize

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
RDLogger.DisableLog("rdApp.*")
torch.set_float32_matmul_precision("medium")


def connectivity(smiles: str) -> str | None:
    molecule = Chem.MolFromSmiles(smiles)
    return Chem.MolToInchiKey(molecule).split("-")[0] if molecule else None


def build_model(task_weights: torch.Tensor):
    message_passing = nn.BondMessagePassing(
        d_h=cfg["message_hidden_dim"],
        depth=cfg["depth"],
        dropout=cfg["mp_dropout"],
        activation=cfg["activation"],
    )
    aggregations = {
        "mean": nn.MeanAggregation,
        "sum": nn.SumAggregation,
        "norm": nn.NormAggregation,
    }
    aggregation = aggregations[cfg["aggregation"]]()
    predictor = nn.RegressionFFN(
        n_tasks=2,
        input_dim=message_passing.output_dim,
        hidden_dim=cfg["ffn_hidden_dim"],
        n_layers=cfg["ffn_num_layers"],
        dropout=cfg["ffn_dropout"],
        criterion=MSE(task_weights=task_weights),
    )
    return models.MPNN(
        message_passing=message_passing,
        agg=aggregation,
        predictor=predictor,
        batch_norm=True,
        warmup_epochs=cfg["warmup_epochs"],
        init_lr=cfg["learning_rate"],
        max_lr=cfg["learning_rate"] * cfg["lr_ratio"],
        final_lr=cfg["learning_rate"] * 0.1,
    )


def loader(smiles: list[str], targets: np.ndarray, shuffle: bool):
    points = [
        chemprop_data.MoleculeDatapoint.from_smi(value, np.asarray(target, dtype=np.float32))
        for value, target in zip(smiles, targets, strict=True)
    ]
    return chemprop_data.build_dataloader(
        chemprop_data.MoleculeDataset(points),
        batch_size=cfg["batch_size"],
        shuffle=shuffle,
    )


def forward(model, smiles: list[str], device: torch.device) -> tuple[np.ndarray, np.ndarray]:
    query_loader = loader(smiles, np.zeros((len(smiles), 2), dtype=np.float32), False)
    predictions = []
    embeddings = []
    model.eval()
    with torch.no_grad():
        for batch in query_loader:
            graph = batch.bmg
            graph.to(device)
            vertex = batch.V_d.to(device) if batch.V_d is not None else None
            extra = batch.X_d.to(device) if batch.X_d is not None else None
            predictions.append(model(graph, vertex, extra).detach().cpu().numpy())
            embeddings.append(model.fingerprint(graph, vertex, extra).detach().cpu().numpy())
    return np.concatenate(predictions), np.concatenate(embeddings)


compounds = pl.read_parquet(inp / cfg["compounds"]).select("compound_id", "standardized_smiles")
challenge_keys = {
    connectivity(smiles) for smiles in compounds["standardized_smiles"] if smiles is not None
}
source = pl.read_parquet(inp / cfg["pxr_single_concentration"])
normalized_cache = {}
records: dict[str, dict] = {}
for smiles, concentration, value in source.select(
    "SMILES", "concentration_M", "log2_fc_estimate"
).iter_rows():
    if smiles not in normalized_cache:
        normalized_cache[smiles] = standardize(smiles)["standardized_smiles"]
    normalized = normalized_cache[smiles]
    key = connectivity(normalized) if normalized else None
    if key is None or key in challenge_keys or value is None:
        continue
    if 8.2e-6 <= concentration <= 8.3e-6:
        head = 0
    elif 3.28e-5 <= concentration <= 3.32e-5:
        head = 1
    else:
        head = None
    if head is None:
        continue
    record = records.setdefault(
        key, {"connectivity": key, "smiles": normalized, "values_0": [], "values_1": []}
    )
    record[f"values_{head}"].append(float(value))

rows = []
for record in records.values():
    values = [record["values_0"], record["values_1"]]
    rows.append(
        {
            "connectivity": record["connectivity"],
            "smiles": record["smiles"],
            "target_0": float(np.mean(values[0])) if values[0] else None,
            "target_1": float(np.mean(values[1])) if values[1] else None,
        }
    )
training = pl.DataFrame(rows).sort("connectivity")
training = training.with_columns(
    pl.struct("smiles", "connectivity")
    .map_elements(
        lambda row: (
            MurckoScaffold.MurckoScaffoldSmiles(smiles=row["smiles"])
            or f"acyclic:{row['connectivity']}"
        ),
        return_dtype=pl.String,
    )
    .alias("scaffold")
)
smiles_train = training["smiles"].to_list()
targets = training.select("target_0", "target_1").to_numpy().astype(np.float32)
means = np.nanmean(targets, axis=0)
stds = np.nanstd(targets, axis=0)
targets_scaled = (targets - means) / stds
task_weights = torch.tensor([1.0, cfg["w_33"]], dtype=torch.float32)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

prediction_sum = np.zeros((compounds.height, 2), dtype=np.float64)
prediction_square_sum = np.zeros_like(prediction_sum)
embedding_sum = np.zeros((compounds.height, cfg["message_hidden_dim"]), dtype=np.float64)
embedding_square_sum = np.zeros_like(embedding_sum)
validation_rows = []
all_smiles = compounds["standardized_smiles"].to_list()
for seed in cfg["seeds"]:
    lightning.seed_everything(seed, workers=True)
    splitter = GroupShuffleSplit(
        n_splits=1, test_size=cfg["validation_fraction"], random_state=seed
    )
    train_index, valid_index = next(
        splitter.split(targets_scaled, groups=training["scaffold"].to_numpy())
    )
    model = build_model(task_weights)
    checkpoint = lightning.callbacks.ModelCheckpoint(
        dirpath=out / f"lightning-seed{seed}",
        filename="best",
        monitor="val_loss",
        mode="min",
        save_top_k=1,
    )
    trainer = lightning.Trainer(
        max_epochs=cfg["max_epochs"],
        accelerator="gpu" if torch.cuda.is_available() else "cpu",
        devices=1,
        enable_progress_bar=False,
        enable_model_summary=False,
        logger=False,
        callbacks=[
            lightning.callbacks.EarlyStopping(
                monitor="val_loss", patience=cfg["patience"], mode="min"
            ),
            checkpoint,
        ],
    )
    trainer.fit(
        model,
        loader(
            [smiles_train[index] for index in train_index],
            targets_scaled[train_index],
            True,
        ),
        loader(
            [smiles_train[index] for index in valid_index],
            targets_scaled[valid_index],
            False,
        ),
    )
    best = torch.load(checkpoint.best_model_path, map_location="cpu", weights_only=False)
    model.load_state_dict(best["state_dict"])
    model.to(device)
    valid_scaled, _ = forward(model, [smiles_train[index] for index in valid_index], device)
    valid_prediction = valid_scaled * stds + means
    for head in range(2):
        observed = targets[valid_index, head]
        mask = np.isfinite(observed)
        validation_rows.append(
            {
                "seed": seed,
                "head": head,
                "n": int(mask.sum()),
                "rmse": float(np.mean((valid_prediction[mask, head] - observed[mask]) ** 2) ** 0.5),
                "spearman": float(
                    spearmanr(valid_prediction[mask, head], observed[mask]).statistic
                ),
                "best_epoch": int(best["epoch"]),
            }
        )
    scaled_prediction, embedding = forward(model, all_smiles, device)
    prediction = scaled_prediction * stds + means
    prediction_sum += prediction
    prediction_square_sum += prediction**2
    embedding_sum += embedding
    embedding_square_sum += embedding**2
    torch.save(
        {
            "state_dict": model.state_dict(),
            "config": cfg,
            "target_means": means.tolist(),
            "target_stds": stds.tolist(),
            "seed": seed,
            "excluded_cyp_connectivities": len(challenge_keys),
        },
        out / f"pxr-chemprop-seed{seed}.pt",
    )
    print(f"seed {seed} complete", flush=True)

n_seeds = len(cfg["seeds"])
prediction_mean = prediction_sum / n_seeds
prediction_std = np.sqrt(np.maximum(prediction_square_sum / n_seeds - prediction_mean**2, 0))
embedding_mean = embedding_sum / n_seeds
embedding_std = np.sqrt(np.maximum(embedding_square_sum / n_seeds - embedding_mean**2, 0))
features = compounds.select("compound_id").with_columns(
    pl.Series("pxr_cp_log2fc_8p25_mean", prediction_mean[:, 0]),
    pl.Series("pxr_cp_log2fc_33_mean", prediction_mean[:, 1]),
    pl.Series("pxr_cp_log2fc_8p25_std", prediction_std[:, 0]),
    pl.Series("pxr_cp_log2fc_33_std", prediction_std[:, 1]),
    pl.Series("pxr_cp_embedding_seed_std", embedding_std.mean(axis=1)),
    *[
        pl.Series(f"pxr_cp_emb_{column:03d}", embedding_mean[:, column])
        for column in range(embedding_mean.shape[1])
    ],
)
features.write_parquet(out / "pxr-chemprop-features.parquet")
pl.DataFrame(validation_rows).write_csv(out / "validation-metrics.csv")
summary = {
    "training_connectivities_after_cyp_exclusion": training.height,
    "excluded_cyp_connectivities": len(challenge_keys),
    "valid_head_0": int(np.isfinite(targets[:, 0]).sum()),
    "valid_head_1": int(np.isfinite(targets[:, 1]).sum()),
    "seeds": cfg["seeds"],
    "device": str(device),
    "validation": validation_rows,
}
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
