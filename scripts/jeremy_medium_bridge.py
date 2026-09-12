"""Bridge upstream medium MP weights to our FT format and verify CPU embeddings."""

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
from finetune_external_multitask import MultitaskModel, encoder_architecture
from rdkit import Chem


def main():
    inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())

    def path(name):
        matches = [p for p in inp.iterdir() if p.name.partition("-")[2] == name]
        assert len(matches) == 1
        return matches[0]

    torch.set_num_threads(2)
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
    original_path = path("chemprop_medium.pt")
    original = torch.load(original_path, map_location="cpu", weights_only=True)
    hyper = original["hyper_parameters"]
    mp_config = dict(hyper["message_passing"])
    assert mp_config.pop("cls") is cnn.BondMessagePassing
    agg_config = dict(hyper["agg"])
    assert agg_config.pop("cls") is cnn.NormAggregation
    assert mp_config["d_h"] == 600 and mp_config["depth"] == 4
    assert agg_config["norm"] == 100 and hyper["batch_norm"] is False
    mp_state = {
        key.removeprefix("message_passing."): tensor
        for key, tensor in original["state_dict"].items()
        if key.startswith("message_passing.")
    }
    bridge = dict(
        state_dict={"mp." + key: tensor for key, tensor in mp_state.items()},
        config=dict(
            architecture="medium",
            width=600,
            depth=4,
            aggregation="norm",
            aggregation_norm=100.0,
            formal_goal_eligible=False,
        ),
        provenance=dict(
            source_repository=cfg["source_repository"],
            source_revision=cfg["source_revision"],
            source_path="CYP_Challenge/checkpoints/chemprop_medium.pt",
            checkpoint_sha256=hashlib.sha256(original_path.read_bytes()).hexdigest(),
            readme_sha256=hashlib.sha256(path("README.md").read_bytes()).hexdigest(),
            role="Upstream public pretrained MP-only bridge; no predictor transferred",
            pretraining_molecule_manifest=None,
            bridge_training_performed=False,
        ),
    )
    torch.save(bridge, out / "best.pt")
    restored = torch.load(out / "best.pt", map_location="cpu", weights_only=True)
    assert all(key.startswith("mp.") for key in restored["state_dict"])
    assert len(restored["state_dict"]) == len(mp_state) == 4
    architecture = encoder_architecture(restored)
    model = MultitaskModel(
        {key.removeprefix("mp."): tensor for key, tensor in restored["state_dict"].items()},
        dropout=0.0,
        n_heads=12,
        architecture=architecture,
    ).eval()
    reference_mp = cnn.BondMessagePassing(**mp_config).eval()
    reference_mp.load_state_dict(mp_state, strict=True)
    reference_aggregation = cnn.NormAggregation(**agg_config).eval()
    compounds = pl.read_parquet(path("compounds.parquet")).sort("compound_id")
    frozen = pl.read_parquet(path("external-embeddings.parquet")).sort("compound_id")
    assert compounds["compound_id"].equals(frozen["compound_id"])
    selection = np.linspace(0, compounds.height - 1, 16, dtype=int)
    selected = compounds[selection.tolist()]
    f = featurizers.SimpleMoleculeMolGraphFeaturizer()
    graph = BatchMolGraph([f(Chem.MolFromSmiles(s)) for s in selected["standardized_smiles"]])
    with torch.inference_mode():
        expected = reference_aggregation(reference_mp(graph), graph.batch).numpy()
        actual = model.embed(graph).numpy()
    original_export = frozen.select(
        sorted(c for c in frozen.columns if c.startswith("emb_"))
    ).to_numpy()[selection]
    assert np.array_equal(actual, expected)
    assert np.allclose(actual, original_export, rtol=1e-5, atol=5e-5)
    selected.select("compound_id").with_columns(
        pl.Series(
            "max_abs_difference_cpu_bridge_vs_original", np.max(np.abs(actual - expected), axis=1)
        ),
        pl.Series(
            "max_abs_difference_cpu_bridge_vs_saved_gpu",
            np.max(np.abs(actual - original_export), axis=1),
        ),
    ).write_parquet(out / "embedding-verification.parquet")
    summary = dict(
        config=cfg,
        bridge_config=restored["config"],
        provenance=restored["provenance"],
        strict_mp_reload=True,
        encoder_tensors=len(mp_state),
        comparison_compounds=len(selection),
        comparison_dimensions=600,
        max_abs_difference_original_cpu=float(np.max(np.abs(actual - expected))),
        max_abs_difference_original_gpu_export=float(np.max(np.abs(actual - original_export))),
        predictor_transferred=False,
        GPU_used=False,
        training_performed=False,
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
