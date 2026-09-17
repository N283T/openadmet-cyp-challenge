"""Interpolate matched joint-encoder weights and audit atomic raw600 transfer."""

import importlib.util
import json
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


class Encoder(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.message_passing = cnn.BondMessagePassing(
            d_h=600, depth=4, dropout=0.0, bias=False
        )
        self.aggregation = cnn.NormAggregation(norm=100.0)

    def forward(self, graph: BatchMolGraph) -> torch.Tensor:
        return self.aggregation(self.message_passing(graph), graph.batch)


torch.set_num_threads(cfg["num_threads"])
torch.backends.cuda.matmul.allow_tf32 = False
device = cfg["device"]
compounds = read("compounds.parquet").sort("compound_id")
assignments = read("parent-assignments.parquet")
measurements = read("measurements.parquet")
spec = importlib.util.spec_from_file_location("official", inp / "custom_scoring_functions.py")
scorer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(scorer)

featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
graphs = []
for smiles in compounds["standardized_smiles"]:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"Invalid standardized SMILES: {smiles}")
    graphs.append(featurizer(molecule))


def graph_batch(indices: np.ndarray) -> BatchMolGraph:
    result = BatchMolGraph([graphs[int(index)] for index in indices])
    result.to(device)
    return result


def encoder_state(name: str) -> dict[str, torch.Tensor]:
    checkpoint = torch.load(inp / name, map_location="cpu", weights_only=False)
    return {
        key.removeprefix("message_passing."): value
        for key, value in checkpoint["state_dict"].items()
        if key.startswith("message_passing.")
    }


def metric(rows: pl.DataFrame, prediction: np.ndarray) -> float:
    return float(
        scorer.rae_soft_threshold_absolute_error(
            rows["value"].to_numpy(),
            prediction,
            y_true_lower=rows["conf_low"].to_numpy(),
            y_true_upper=rows["conf_high"].to_numpy(),
        )
    )


results = []
for pair in cfg["pairs"]:
    base = encoder_state(pair["base"])
    focused = encoder_state(pair["focused"])
    if base.keys() != focused.keys():
        raise ValueError(f"Encoder keys differ for {pair['name']}")
    table = (
        measurements.filter(pl.col("endpoint") == pair["endpoint"])
        .join(assignments, on="compound_id", validate="1:1")
        .sort("compound_id")
    )
    for alpha in cfg["alphas"]:
        method = f"{pair['name']}_alpha{int(round(alpha * 100)):03d}"
        state = {
            key: (1.0 - alpha) * base[key] + alpha * focused[key] for key in base
        }
        encoder = Encoder().to(device)
        encoder.message_passing.load_state_dict(state)
        encoder.eval()
        chunks = []
        with torch.inference_mode():
            for start in range(0, compounds.height, cfg["batch_size"]):
                indices = np.arange(start, min(start + cfg["batch_size"], compounds.height))
                chunks.append(encoder(graph_batch(indices)).cpu().numpy())
        values = np.concatenate(chunks).astype(np.float32)
        embeddings = pl.DataFrame(values, schema=[f"emb_{i:03d}" for i in range(600)])
        embeddings.insert_column(0, compounds["compound_id"])
        trial_dir = out / method
        trial_dir.mkdir()
        embeddings.write_parquet(trial_dir / "external-embeddings.parquet")
        torch.save(
            {
                "message_passing_state_dict": state,
                "base": pair["base"],
                "focused": pair["focused"],
                "alpha": alpha,
            },
            trial_dir / "interpolated-encoder.pt",
        )
        joined = table.join(embeddings, on="compound_id", validate="1:1")
        columns = [f"emb_{i:03d}" for i in range(600)]
        x = joined.select(columns).to_numpy()
        y = joined["value"].to_numpy()
        predictions = []
        for fold in range(5):
            valid = (joined["fold"].to_numpy() == fold) & ~joined["is_anchor"].to_numpy()
            fitted = lgb.LGBMRegressor(
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
            fitted.fit(x[~valid], y[~valid])
            selected = joined.filter(pl.Series(valid))
            predictions.append(
                selected.select(
                    "compound_id",
                    "endpoint",
                    "fold",
                    "is_neighbor",
                    "conf_low",
                    "conf_high",
                    "value",
                ).with_columns(
                    pl.lit(method).alias("method"),
                    pl.Series("prediction", fitted.predict(x[valid]).astype(float)),
                )
            )
        oof = pl.concat(predictions)
        oof.write_parquet(trial_dir / "external-oof.parquet")
        result = {
            "pair": pair["name"],
            "endpoint": pair["endpoint"],
            "alpha": alpha,
            "method": method,
            "ST_RAE": metric(oof, oof["prediction"].to_numpy()),
        }
        results.append(result)
        pl.DataFrame(results).write_csv(out / "metrics-progress.csv")
        print(json.dumps(result), flush=True)

pl.DataFrame(results).write_csv(out / "metrics.csv")
(out / "summary.json").write_text(json.dumps({"config": cfg, "results": results}, indent=2) + "\n")
