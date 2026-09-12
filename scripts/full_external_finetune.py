"""Refit the adopted Jeremy-medium encoder on all available challenge labels.

No validation subset is reserved. Epoch count is the median one-based selected
epoch from the declared, already-completed OOF fine-tuning run.
"""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
from chemprop import featurizers
from chemprop.data import BatchMolGraph
from finetune_external_multitask import (
    MultitaskModel,
    encoder_architecture,
    input_path,
    masked_loss,
)
from rdkit import Chem


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    previous = json.loads(input_path(inp, "summary.json").read_text())
    recipe = previous["config"]
    selected_epochs = [fold["best_epoch"] + 1 for fold in previous["folds"]]
    epochs = int(np.median(selected_epochs))
    assert len(selected_epochs) == 5 and epochs > 0
    assert recipe["training_loss"] == "mse"
    cfg.update(
        epochs=epochs,
        selected_one_based_epochs=selected_epochs,
        epoch_policy="Median one-based OOF-selected epoch; fixed full-data training",
        recipe=recipe,
    )
    (out / "effective-config.json").write_text(json.dumps(cfg, indent=2))
    compounds = pl.read_parquet(input_path(inp, "compounds.parquet")).sort("compound_id")
    assignments = pl.read_parquet(input_path(inp, "parent-assignments.parquet"))
    labels = pl.read_parquet(input_path(inp, "measurements.parquet")).filter(
        pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition")
    )
    auxiliary = pl.read_parquet(input_path(inp, "targets.parquet"))
    endpoints, aux_heads = sorted(labels["endpoint"].unique()), sorted(auxiliary["head"].unique())
    assert len(endpoints) == len(aux_heads) == 4
    heads = endpoints + aux_heads
    targets = pl.concat(
        [
            labels.select("compound_id", pl.col("endpoint").alias("head"), "value"),
            auxiliary.select("compound_id", "head", "value"),
        ]
    )
    assert targets.select("compound_id", "head").unique().height == targets.height
    assert set(targets["compound_id"]).issubset(set(assignments["compound_id"]))
    assert compounds.height == compounds["compound_id"].n_unique() == 6896
    assert assignments.height == assignments["compound_id"].n_unique() == 6146
    assert labels.height == 6525 and auxiliary.height == 17504
    h = compounds.select("compound_id").join(
        targets.pivot(on="head", index="compound_id", values="value"),
        on="compound_id",
        how="left",
        validate="1:1",
        maintain_order="left",
    )
    y = h.select(heads).to_numpy().astype(np.float32)
    training = np.flatnonzero(np.isfinite(y).any(1))
    non_test = np.isin(compounds["compound_id"].to_numpy(), assignments["compound_id"].to_numpy())
    assert non_test.sum() == 6146 and not np.isfinite(y[~non_test]).any()
    assert set(h["compound_id"][training]).issubset(set(assignments["compound_id"]))
    means, stds = np.nanmean(y[training], axis=0), np.nanstd(y[training], axis=0)
    assert np.isfinite(means).all() and np.isfinite(stds).all()
    stds = np.where(stds < 1e-6, 1.0, stds)
    z = (y - means) / stds
    counts = dict(zip(heads, np.isfinite(y[training]).sum(0).tolist(), strict=True))
    roles = np.where(non_test, "no_available_labels", "test_inference_only")
    roles[np.isfinite(y).any(1)] = "full_training"
    compounds.select("compound_id").with_columns(
        pl.Series("role", roles),
        pl.Series("available_heads", np.isfinite(y).sum(1)),
        pl.lit("full_training").alias("split_id"),
    ).write_parquet(out / "membership.parquet")
    targets.with_columns(pl.lit("full_training").alias("split_id")).write_parquet(
        out / "training-labels.parquet"
    )
    torch.set_num_threads(recipe["num_threads"])
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.manual_seed(recipe["seed"])
    rng = np.random.default_rng(recipe["seed"])
    checkpoint = torch.load(input_path(inp, "best.pt"), map_location="cpu", weights_only=False)
    encoder_state = {
        key.removeprefix("mp."): value
        for key, value in checkpoint["state_dict"].items()
        if key.startswith("mp.")
    }
    architecture = encoder_architecture(checkpoint)
    assert architecture == previous["encoder_architecture"]
    device = recipe["device"]
    model = MultitaskModel(encoder_state, recipe["dropout"], len(heads), architecture).to(device)
    optimizer = torch.optim.AdamW(
        [
            {"params": model.mp.parameters(), "lr": recipe["encoder_lr"]},
            {"params": model.head.parameters(), "lr": recipe["head_lr"]},
        ],
        weight_decay=recipe["weight_decay"],
    )
    weights = torch.tensor(
        [recipe["pIC50_weight"]] * 4 + [recipe["log2fc_weight"]] * 4,
        dtype=torch.float32,
        device=device,
    )
    featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
    graphs = []
    for smiles in compounds["standardized_smiles"]:
        molecule = Chem.MolFromSmiles(smiles)
        assert molecule is not None
        graphs.append(featurizer(molecule))

    def graph(indices):
        batch = BatchMolGraph([graphs[int(i)] for i in indices])
        batch.to(device)
        return batch

    for epoch in range(epochs):
        model.train()
        total, count = 0.0, 0
        order = rng.permutation(training)
        for start in range(0, len(order), recipe["batch_size"]):
            batch = order[start : start + recipe["batch_size"]]
            optimizer.zero_grad(set_to_none=True)
            loss = masked_loss(
                model(graph(batch)), torch.as_tensor(z[batch], device=device), weights, "mse"
            )
            assert torch.isfinite(loss)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0, error_if_nonfinite=True)
            optimizer.step()
            total += float(loss.detach()) * len(batch)
            count += len(batch)
        record = dict(epoch=epoch, train_mse=total / count, training_compounds=count)
        assert count == len(training)
        with (out / "history.jsonl").open("a") as stream:
            stream.write(json.dumps(record) + "\n")
        print(json.dumps(record), flush=True)
    state = dict(
        state_dict=model.state_dict(),
        optimizer=optimizer.state_dict(),
        epoch=epochs - 1,
        heads=heads,
        means=means.tolist(),
        stds=stds.tolist(),
        config=cfg,
        encoder_architecture=architecture,
        split_id="full_training",
    )
    torch.save(state, out / "full-model.pt")
    model.eval()
    embeddings, predictions = [], []
    with torch.inference_mode():
        for start in range(0, compounds.height, recipe["inference_batch_size"]):
            batch = np.arange(start, min(start + recipe["inference_batch_size"], compounds.height))
            emb = model.embed(graph(batch))
            embeddings.append(emb.cpu().numpy())
            predictions.append(model.head(emb).cpu().numpy() * stds + means)
    embeddings, predictions = np.concatenate(embeddings), np.concatenate(predictions)
    assert embeddings.shape == (6896, 600) and predictions.shape == (6896, 8)
    assert np.isfinite(embeddings).all() and np.isfinite(predictions).all()
    for filename, matrix, names in [
        ("external-embeddings.parquet", embeddings, [f"emb_{i:03d}" for i in range(600)]),
        ("direct-predictions.parquet", predictions, heads),
    ]:
        pl.DataFrame(matrix, schema=names).insert_column(0, compounds["compound_id"]).with_columns(
            pl.lit("full_training").alias("split_id"), pl.lit(-1).alias("fold")
        ).write_parquet(out / filename)
    summary = dict(
        epochs=epochs,
        selected_one_based_epochs=selected_epochs,
        training_compounds=len(training),
        non_test_compounds=int(non_test.sum()),
        test_inference_only=int((~non_test).sum()),
        training_labels=counts,
        test_labels_used=0,
        reserved_validation_compounds=0,
        dimensions=600,
        exported_compounds=compounds.height,
        encoder_architecture=architecture,
        source_finetune_run=cfg["source_finetune_run"],
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
