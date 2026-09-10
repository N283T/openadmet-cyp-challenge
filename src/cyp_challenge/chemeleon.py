"""Frozen CheMeleon raw message-passing fingerprints with explicit molecule keys."""

import importlib.metadata
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
from chemprop import featurizers, nn
from chemprop.data import BatchMolGraph
from chemprop.models import MPNN
from rdkit import Chem

from .embedding_validation import validate_embeddings
from .storage import digest, write_json


class Fingerprint:
    def __init__(self, checkpoint, device):
        saved = torch.load(checkpoint, map_location="cpu", weights_only=True)
        self.mp = nn.BondMessagePassing(**saved["hyper_parameters"])
        self.mp.load_state_dict(saved["state_dict"], strict=True)
        self.mp.to(device).eval().requires_grad_(False)
        self.aggregation = nn.MeanAggregation().to(device).eval()
        self.featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
        self.device = device
        self.hyperparameters = saved["hyper_parameters"]

    def graph(self, smiles):
        mols = [Chem.MolFromSmiles(s) for s in smiles]
        if any(m is None for m in mols):
            raise ValueError("Invalid standardized SMILES")
        graph = BatchMolGraph([self.featurizer(m) for m in mols])
        graph.to(self.device)
        return graph

    def __call__(self, smiles):
        graph = self.graph(smiles)
        with torch.inference_mode():
            values = self.aggregation(self.mp(graph), graph.batch).cpu().numpy()
        if values.shape != (len(smiles), self.mp.output_dim) or not np.isfinite(values).all():
            raise ValueError("Invalid CheMeleon output shape or nonfinite values")
        return values.astype(np.float32, copy=False)


def main():
    inputs, output = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
    torch.manual_seed(config["seed"])
    np.random.seed(config["seed"])
    torch.set_num_threads(config["threads"])
    torch.use_deterministic_algorithms(True)
    torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False
    checkpoint = next(inputs.glob("*-chemeleon_mp.pt"))
    provenance = json.loads(next(inputs.glob("*-provenance.json")).read_text())
    assert digest(checkpoint) == provenance["checkpoint_sha256"]
    compounds = pl.read_parquet(next(inputs.glob("*-compounds.parquet"))).sort("compound_id")
    assert compounds["compound_id"].n_unique() == compounds.height
    encoder = Fingerprint(checkpoint, config["device"])
    assert encoder.mp.output_dim == config["dimension"]
    runtime = {
        "config": config,
        "checkpoint_sha256": digest(checkpoint),
        "versions": {
            name: importlib.metadata.version(name)
            for name in ["torch", "chemprop", "rdkit", "numpy", "polars"]
        },
        "cuda_runtime": torch.version.cuda,
        "device": torch.cuda.get_device_name() if config["device"] == "cuda" else "cpu",
        "hyperparameters": encoder.hyperparameters,
    }
    write_json(output / "runtime.json", runtime)
    smiles = compounds["standardized_smiles"].to_list()
    frames = []
    names = [f"chemeleon_{i:04d}" for i in range(config["dimension"])]
    for start in range(0, compounds.height, config["batch_size"]):
        batch = compounds.slice(start, config["batch_size"])
        values = encoder(batch["standardized_smiles"].to_list())
        frame = pl.DataFrame(values, schema=names).insert_column(0, batch["compound_id"])
        frame.write_parquet(output / f"chemeleon-{start:05d}.parquet")
        frames.append(frame)
        if start % 1024 == 0:
            print(f"{start + batch.height}/{compounds.height}", flush=True)
    features = pl.concat(frames)
    validate_embeddings(features, compounds["compound_id"].to_list(), names)
    features.write_parquet(output / "embeddings.parquet")
    probe = smiles[:8]
    first = features.head(8).select(names).to_numpy()
    repeated = encoder(probe)
    singleton = np.concatenate([encoder([s]) for s in probe])
    # Match the official MPNN.fingerprint path, which never uses predictor weights.
    reference = (
        MPNN(encoder.mp, encoder.aggregation, nn.RegressionFFN(input_dim=encoder.mp.output_dim))
        .to(config["device"])
        .eval()
    )
    with torch.inference_mode():
        official_path = reference.fingerprint(encoder.graph(probe)).cpu().numpy()
    tolerance = config["validation_tolerance"]
    for other in [repeated, singleton, official_path]:
        np.testing.assert_allclose(first, other, **tolerance)
    write_json(
        output / "summary.json",
        {
            **runtime,
            "rows": features.height,
            "dimension": len(names),
            "columns": names,
            "output_sha256": digest(output / "embeddings.parquet"),
            "validation": {
                "n": len(probe),
                "repeat_max_abs": float(np.abs(first - repeated).max()),
                "singleton_max_abs": float(np.abs(first - singleton).max()),
                "official_path_max_abs": float(np.abs(first - official_path).max()),
                "finite": True,
                "unique_complete_compound_ids": True,
            },
            "training": "none; shared frozen features, no CYP labels consumed",
        },
    )
    print("Completed", features.shape, flush=True)


if __name__ == "__main__":
    main()
