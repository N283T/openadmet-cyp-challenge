"""CPU-only twelve-head preparation and aggregation compatibility checks."""

import json
import os
from pathlib import Path

import finetune_external_multitask as ft
import torch
from chemprop import featurizers
from chemprop.data import BatchMolGraph
from rdkit import Chem

ft.main()
inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
checkpoint = torch.load(ft.input_path(inp, "best.pt"), map_location="cpu", weights_only=False)
state = {
    k.removeprefix("mp."): v for k, v in checkpoint["state_dict"].items() if k.startswith("mp.")
}
summary = json.loads((out / "summary.json").read_text())
assert len(summary["heads"]) == 12 and not summary["training_executed"]
architecture = ft.encoder_architecture(checkpoint)
model = ft.MultitaskModel(state, 0.1, n_heads=12, architecture=architecture).eval()
default = ft.MultitaskModel(state, 0.1, n_heads=12).eval()
f = featurizers.SimpleMoleculeMolGraphFeaturizer()
g = BatchMolGraph([f(Chem.MolFromSmiles(s)) for s in ["C", "CCO"]])
with torch.inference_mode():
    mh = model.mp(g)
    expected = torch.zeros((2, 600)).index_add_(0, g.batch, mh) / 100
    assert torch.allclose(model.embed(g), expected, atol=1e-7, rtol=1e-6)
    dh = default.mp(g)
    expected_mean = torch.zeros((2, 600)).index_add_(0, g.batch, dh) / torch.bincount(
        g.batch
    ).reshape(-1, 1)
    assert torch.allclose(default.embed(g), expected_mean, atol=1e-7, rtol=1e-6)
    assert model(g).shape == default(g).shape == (2, 12)
assert model.mp.depth == 4 and default.mp.depth == 5
assert model.head[0].out_features == default.head[0].out_features == 256
(out / "model-smoke.json").write_text(
    json.dumps(
        dict(
            medium=architecture,
            default=ft.encoder_architecture({}),
            strict_mp_transfer=True,
            medium_sum_norm100=True,
            default_mean=True,
            output_heads=12,
            challenge_ffn_hidden_width=256,
            training_executed=False,
        ),
        indent=2,
    )
)
print("Twelve-head medium/default CPU model smoke passed", flush=True)
