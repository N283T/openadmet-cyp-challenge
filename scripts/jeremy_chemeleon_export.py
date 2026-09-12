"""Frozen extraction from Jeremy's upstream public chemprop_chemeleon checkpoint."""

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
from chemprop import featurizers
from chemprop import nn as cnn
from chemprop.data import BatchMolGraph
from chemprop.nn.transforms import GraphTransform, UnscaleTransform
from rdkit import Chem


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

    def path(name):
        matches = [p for p in inp.iterdir() if p.name.partition("-")[2] == name]
        assert len(matches) == 1
        return matches[0]

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
    checkpoint_path = path("chemprop_chemeleon.pt")
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
    hyper = checkpoint["hyper_parameters"]
    mp_config = dict(hyper["message_passing"])
    assert mp_config.pop("cls") is cnn.BondMessagePassing
    agg_config = dict(hyper["agg"])
    assert agg_config.pop("cls") is cnn.MeanAggregation
    assert mp_config["d_h"] == 2048 and mp_config["depth"] == 6
    assert hyper["batch_norm"] is False
    assert hyper["predictor"]["n_tasks"] == 47
    torch.set_num_threads(2)
    torch.backends.cuda.matmul.allow_tf32 = False
    device = cfg.get("device", "cuda")
    mp = cnn.BondMessagePassing(**mp_config).to(device).eval()
    mp.load_state_dict(
        {
            key.removeprefix("message_passing."): value
            for key, value in checkpoint["state_dict"].items()
            if key.startswith("message_passing.")
        },
        strict=True,
    )
    agg = cnn.MeanAggregation(**agg_config).to(device).eval()
    compounds = pl.read_parquet(path("compounds.parquet")).sort("compound_id")
    assert compounds["compound_id"].n_unique() == compounds.height == 6896
    featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
    values = []
    with torch.inference_mode():
        for start in range(0, compounds.height, 64):
            graphs = []
            for smiles in compounds.slice(start, 64)["standardized_smiles"]:
                mol = Chem.MolFromSmiles(smiles)
                assert mol is not None
                graphs.append(featurizer(mol))
            batch = BatchMolGraph(graphs)
            batch.to(device)
            values.append(agg(mp(batch), batch.batch).cpu().numpy())
    values = np.concatenate(values)
    assert values.shape == (6896, 2048) and np.isfinite(values).all()
    pl.DataFrame(values, schema=[f"emb_{i:03d}" for i in range(2048)]).insert_column(
        0, compounds["compound_id"]
    ).write_parquet(out / "external-embeddings.parquet")
    provenance = dict(
        config=cfg,
        checkpoint_sha256=hashlib.sha256(checkpoint_path.read_bytes()).hexdigest(),
        checkpoint_role="Upstream sibling-project pretraining checkpoint per README and script03",
        readme_sha256=hashlib.sha256(path("README.md").read_bytes()).hexdigest(),
        extraction_source_sha256=hashlib.sha256(
            path("03_embed_frozen_encoders.py").read_bytes()
        ).hexdigest(),
        dimensions=2048,
        depth=6,
        aggregation="MeanAggregation",
        graph_transform=repr(mp_config.get("graph_transform")),
        pretraining_output_heads=47,
        pretraining_task_names=None,
        pretraining_molecule_manifest=None,
        formal_goal_eligible=False,
        rows=compounds.height,
        checkpoint_is_challenge_17_head_finetune=False,
        challenge_labels_used_during_extraction=False,
        upstream_overlap_status="Unknown: upstream training molecule manifest not provided",
        modification="None: strict encoder weights, checkpoint architecture and frozen aggregation",
    )
    (out / "summary.json").write_text(json.dumps(provenance, indent=2))
    print(json.dumps(provenance), flush=True)


if __name__ == "__main__":
    main()
