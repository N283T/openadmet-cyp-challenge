"""Frozen upstream ADME predictor features; no challenge labels or fitting."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
from chemprop import nn as cnn
from chemprop.nn.transforms import GraphTransform, UnscaleTransform

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def path(name):
    matches = [p for p in inp.iterdir() if p.name.partition("-")[2] == name]
    assert len(matches) == 1, name
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
torch.set_num_threads(3)
checkpoint = torch.load(path(cfg["checkpoint_file"]), map_location="cpu", weights_only=True)
hyper = dict(checkpoint["hyper_parameters"]["predictor"])
assert hyper.pop("cls") is cnn.RegressionFFN
assert hyper["n_tasks"] == 47 and hyper["n_layers"] == 3
predictor = cnn.RegressionFFN(**hyper).eval()
predictor.load_state_dict(
    {
        k.removeprefix("predictor."): v
        for k, v in checkpoint["state_dict"].items()
        if k.startswith("predictor.")
    },
    strict=True,
)
frame = pl.read_parquet(path("external-embeddings.parquet")).sort("compound_id")
columns = sorted(c for c in frame.columns if c.startswith("emb_"))
x = torch.as_tensor(frame.select(columns).to_numpy(), dtype=torch.float32)
assert frame["compound_id"].n_unique() == frame.height == 6896
with torch.inference_mode():
    hidden = predictor.ffn[-1][:-1](predictor.ffn[:-1](x))
    standardized = predictor.ffn[-1][-1](hidden)
    torch.testing.assert_close(standardized, predictor.ffn(x), atol=1e-6, rtol=1e-6)
    scores = predictor(x)
assert hidden.shape[0] == scores.shape[0] == 6896 and scores.shape[1] == 47
for values, prefix, filename in [
    (hidden, "emb_", "last-hidden.parquet"),
    (scores, "score_", "external-scores.parquet"),
]:
    array = values.numpy()
    assert np.isfinite(array).all()
    pl.DataFrame(array, schema=[f"{prefix}{i:03d}" for i in range(array.shape[1])]).insert_column(
        0, frame["compound_id"]
    ).write_parquet(out / filename)
(out / "summary.json").write_text(
    json.dumps(
        dict(
            config=cfg,
            rows=frame.height,
            hidden_dimensions=hidden.shape[1],
            score_dimensions=47,
            challenge_labels_used=False,
            fitting=False,
            task_names=None,
            source_role="Original upstream 47-task predictor, strictly loaded and evaluated",
            feature_semantics="Final activated hidden representation and predictor outputs",
            upstream_manifest="Not supplied; no invented task names",
        ),
        indent=2,
    )
)
