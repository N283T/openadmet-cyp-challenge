"""Independently check full-fit membership, label coverage, epochs and saved export."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
from chemprop import featurizers
from chemprop.data import BatchMolGraph
from finetune_external_multitask import MultitaskModel, input_path
from rdkit import Chem


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    compounds = pl.read_parquet(input_path(inp, "compounds.parquet")).sort("compound_id")
    members = pl.read_parquet(input_path(inp, "membership.parquet")).sort("compound_id")
    embeddings = pl.read_parquet(input_path(inp, "external-embeddings.parquet")).sort("compound_id")
    predictions = pl.read_parquet(input_path(inp, "direct-predictions.parquet")).sort("compound_id")
    saved_labels = pl.read_parquet(input_path(inp, "training-labels.parquet"))
    pic50 = (
        pl.read_parquet(input_path(inp, "measurements.parquet"))
        .filter(pl.col("endpoint").str.ends_with("_pIC50_direct_inhibition"))
        .select("compound_id", pl.col("endpoint").alias("head"), "value")
    )
    auxiliary = pl.read_parquet(input_path(inp, "targets.parquet")).select(
        "compound_id", "head", "value"
    )
    expected_labels = pl.concat([pic50, auxiliary]).sort("compound_id", "head")
    assert (
        saved_labels.select(expected_labels.columns)
        .sort("compound_id", "head")
        .equals(expected_labels)
    )
    assert expected_labels.height == 24029
    for table in [members, embeddings, predictions]:
        assert table["compound_id"].n_unique() == 6896
        assert table["compound_id"].equals(compounds["compound_id"])
    assert members.filter(pl.col("role") == "test_inference_only").height == 750
    training_ids = set(members.filter(pl.col("role") == "full_training")["compound_id"])
    assert training_ids == set(expected_labels["compound_id"])
    assert len(training_ids) == 4906
    checkpoint = torch.load(
        input_path(inp, "full-model.pt"), map_location="cpu", weights_only=False
    )
    selected = json.loads(input_path(inp, "summary.json").read_text())
    expected_epochs = int(np.median([fold["best_epoch"] + 1 for fold in selected["folds"]]))
    assert checkpoint["epoch"] + 1 == expected_epochs == 94
    heads = checkpoint["heads"]
    y = (
        compounds.select("compound_id")
        .join(
            expected_labels.pivot(on="head", index="compound_id", values="value"),
            on="compound_id",
            how="left",
            validate="1:1",
            maintain_order="left",
        )
        .select(heads)
        .to_numpy()
        .astype(np.float32)
    )
    fitted_y = y[np.isfinite(y).any(1)]
    assert np.array_equal(np.nanmean(fitted_y, axis=0), checkpoint["means"])
    assert np.array_equal(np.nanstd(fitted_y, axis=0), checkpoint["stds"])
    x = embeddings.select([f"emb_{i:03d}" for i in range(600)]).to_numpy()
    pred = predictions.select(heads).to_numpy()
    assert np.isfinite(x).all() and np.isfinite(pred).all()
    model = MultitaskModel(
        {
            k.removeprefix("mp."): v
            for k, v in checkpoint["state_dict"].items()
            if k.startswith("mp.")
        },
        selected["config"]["dropout"],
        len(heads),
        checkpoint["encoder_architecture"],
    )
    model.load_state_dict(checkpoint["state_dict"], strict=True)
    model.eval()
    torch.set_num_threads(4)
    indices = np.linspace(0, 6895, 16, dtype=int)
    featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
    graph = BatchMolGraph(
        [featurizer(Chem.MolFromSmiles(compounds["standardized_smiles"][int(i)])) for i in indices]
    )
    with torch.inference_mode():
        actual_x = model.embed(graph).numpy()
        actual_pred = model.head(torch.tensor(actual_x)).numpy()
        actual_pred = actual_pred * checkpoint["stds"] + checkpoint["means"]
    assert np.allclose(actual_x, x[indices], rtol=1e-4, atol=1e-5)
    assert np.allclose(actual_pred, pred[indices], rtol=1e-4, atol=1e-5)
    result = dict(
        passed=True,
        epochs=expected_epochs,
        training_compounds=len(training_ids),
        training_labels=expected_labels.height,
        test_inference_only=750,
        label_table_exact=True,
        all_training_scalers_exact=True,
        model_export_probe_rows=len(indices),
        embedding_max_abs_difference=float(np.max(np.abs(actual_x - x[indices]))),
        prediction_max_abs_difference=float(np.max(np.abs(actual_pred - pred[indices]))),
    )
    (out / "verification.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
