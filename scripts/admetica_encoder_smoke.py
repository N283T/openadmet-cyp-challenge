"""Load only reviewed checkpoint types and verify reusable CYP encoder weights."""

import json
import os
from pathlib import Path

import numpy as np
import torch
from chemprop import featurizers, nn
from chemprop.data import BatchMolGraph
from rdkit import Chem

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
allowed = [
    nn.MeanAggregation,
    nn.BondMessagePassing,
    nn.BinaryClassificationFFN,
    nn.Activation,
    (dict, "lightning.fabric.utilities.data.AttributeDict"),
    (nn.BinaryAccuracy, "chemprop.nn.metrics.BinaryAccuracyMetric"),
    (nn.BinaryAUROC, "chemprop.nn.metrics.BinaryAUROCMetric"),
    (np._core.multiarray.scalar, "numpy.core.multiarray.scalar"),
    (np._core.multiarray._reconstruct, "numpy.core.multiarray._reconstruct"),
    np.ndarray,
    np.dtype,
    np.dtypes.Float64DType,
]
results = []
for p in sorted(inp.glob("*.ckpt")):
    with torch.serialization.safe_globals(allowed):
        ckpt = torch.load(p, map_location="cpu", weights_only=True)
    params = ckpt["hyper_parameters"]["message_passing"]
    params = {k: v for k, v in params.items() if k != "cls"}
    mp = nn.BondMessagePassing(**params).eval()
    state = {
        k.removeprefix("message_passing."): v
        for k, v in ckpt["state_dict"].items()
        if k.startswith("message_passing.")
    }
    mp.load_state_dict(state, strict=True)
    smiles = ["CCO", "c1ccccc1", "Cn1c(=O)c2c(ncn2C)n(C)c1=O"]
    f = featurizers.SimpleMoleculeMolGraphFeaturizer()
    graph = BatchMolGraph([f(Chem.MolFromSmiles(s)) for s in smiles])
    with torch.inference_mode():
        emb = nn.MeanAggregation()(mp(graph), graph.batch).numpy()
    assert emb.shape == (3, 300) and np.isfinite(emb).all()
    torch.save(dict(hyper_parameters=params, state_dict=state), out / (p.stem + "-encoder.pt"))
    results.append(
        dict(
            file=p.name,
            shape=list(emb.shape),
            finite=True,
            strict_state_match=True,
            parameters={k: str(v) for k, v in params.items()},
        )
    )
(out / "summary.json").write_text(json.dumps(results, indent=2))
print(json.dumps(results, indent=2))
