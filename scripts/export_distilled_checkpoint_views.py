"""Export raw, predictor-hidden, and predictor-score views from owned checkpoints."""

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

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


class JointModel(torch.nn.Module):
    def __init__(self, config: dict, n_tasks: int):
        super().__init__()
        self.message_passing = cnn.BondMessagePassing(
            d_h=config["width"], depth=config["depth"], dropout=0.0, bias=False
        )
        self.aggregation = cnn.NormAggregation(norm=config["aggregation_norm"])
        layers = []
        for _ in range(config["ffn_hidden_layers"]):
            layers.extend(
                [
                    torch.nn.Linear(config["width"], config["width"]),
                    torch.nn.ReLU(),
                    torch.nn.Dropout(config["dropout"]),
                ]
            )
        layers.append(torch.nn.Linear(config["width"], n_tasks))
        self.predictor = torch.nn.Sequential(*layers)

    def embed(self, graph: BatchMolGraph) -> torch.Tensor:
        return self.aggregation(self.message_passing(graph), graph.batch)


compounds = pl.read_parquet(inp / "compounds.parquet").sort("compound_id")
assert compounds["compound_id"].n_unique() == compounds.height
featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
graphs = []
for smiles in compounds["standardized_smiles"]:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Invalid standardized SMILES: {smiles}")
    graphs.append(featurizer(molecule))

torch.set_num_threads(cfg["num_threads"])
device = cfg["device"]
summary = []
for item in cfg["checkpoints"]:
    checkpoint = torch.load(inp / item["file"], map_location="cpu", weights_only=False)
    model_config = checkpoint["config"]
    tasks = checkpoint["tasks"]
    model = JointModel(model_config, len(tasks))
    model.load_state_dict(checkpoint["state_dict"], strict=True)
    model.to(device).eval()
    raw_chunks, hidden_chunks, score_chunks = [], [], []
    with torch.inference_mode():
        for start in range(0, len(graphs), cfg["batch_size"]):
            batch = BatchMolGraph(graphs[start : start + cfg["batch_size"]])
            batch.to(device)
            raw = model.embed(batch)
            hidden = model.predictor[:-1](raw)
            scores = model.predictor[-1](hidden)
            raw_chunks.append(raw.cpu().numpy())
            hidden_chunks.append(hidden.cpu().numpy())
            score_chunks.append(scores.cpu().numpy())
    raw = np.concatenate(raw_chunks).astype(np.float32)
    hidden = np.concatenate(hidden_chunks).astype(np.float32)
    scores = np.concatenate(score_chunks).astype(np.float32)
    for values in (raw, hidden, scores):
        if not np.isfinite(values).all():
            raise ValueError(f"Non-finite values in {item['name']}")

    def frame(values: np.ndarray, prefix: str) -> pl.DataFrame:
        result = pl.DataFrame(
            values, schema=[f"{prefix}{index:03d}" for index in range(values.shape[1])]
        )
        result.insert_column(0, compounds["compound_id"])
        return result

    raw_frame = frame(raw, "raw_")
    raw_frame.write_parquet(out / f"{item['name']}-raw.parquet")
    exported = ["raw"]
    if item["predictor_trained"]:
        hidden_frame = frame(hidden, "hidden_")
        score_frame = frame(scores, "score_")
        hidden_frame.write_parquet(out / f"{item['name']}-hidden.parquet")
        score_frame.write_parquet(out / f"{item['name']}-scores.parquet")
        raw_frame.join(score_frame, on="compound_id", validate="1:1").write_parquet(
            out / f"{item['name']}-raw-scores.parquet"
        )
        exported.extend(["hidden", "scores", "raw_scores"])
    summary.append(
        {
            "name": item["name"],
            "tasks": len(tasks),
            "task_ids": [task["assay_id"] for task in tasks],
            "predictor_trained": item["predictor_trained"],
            "views": exported,
            "raw_dimensions": raw.shape[1],
            "hidden_dimensions": hidden.shape[1],
            "score_dimensions": scores.shape[1],
        }
    )
    print(json.dumps(summary[-1]), flush=True)
    del model
    torch.cuda.empty_cache()

(out / "summary.json").write_text(json.dumps({"checkpoints": summary}, indent=2))
