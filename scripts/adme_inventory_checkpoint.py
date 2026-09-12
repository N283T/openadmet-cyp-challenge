"""Inspect checkpoint metadata without executing untrusted pickle globals or training."""

import json
import os
import pickletools
import zipfile
from pathlib import Path

import torch
from chemprop import nn as cnn
from chemprop.nn.transforms import GraphTransform, UnscaleTransform

torch.set_num_threads(1)
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


def summarize(x, depth=0):
    if isinstance(x, torch.Tensor):
        return dict(shape=list(x.shape), dtype=str(x.dtype))
    if isinstance(x, dict):
        return {str(k): summarize(v, depth + 1) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [summarize(v, depth + 1) for v in x]
    if isinstance(x, (str, int, float, bool)) or x is None:
        return x
    return str(type(x))


results = {}
for p in Path(os.environ["CYP_INPUT_DIR"]).iterdir():
    if p.suffix == ".pt":
        try:
            obj = torch.load(p, map_location="cpu", weights_only=True)
            results[p.name] = summarize(obj)
        except Exception as e:
            globals_found = []
            strings = []
            with zipfile.ZipFile(p) as z:
                f = next(n for n in z.namelist() if n.endswith("data.pkl"))
                for op, arg, _ in pickletools.genops(z.read(f)):
                    if op.name == "GLOBAL":
                        globals_found.append(arg)
                    if op.name in ["BINUNICODE", "SHORT_BINUNICODE"] and isinstance(arg, str):
                        strings.append(arg)
            results[p.name] = dict(safe_load_error=str(e), globals=globals_found, strings=strings)
    elif p.suffix == ".json":
        results[p.name] = json.loads(p.read_text())
Path(os.environ["CYP_OUTPUT_DIR"], "metadata.json").write_text(json.dumps(results, indent=2))
print(
    json.dumps(
        {k: list(v) if isinstance(v, dict) else type(v).__name__ for k, v in results.items()}
    )
)
