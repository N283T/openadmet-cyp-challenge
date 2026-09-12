"""Export trained broad-panel task scores from saved frozen embeddings and head."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
checkpoint = torch.load(next(inp.glob("*-best.pt")), map_location="cpu", weights_only=False)
emb = pl.read_parquet(next(inp.glob("*-external-embeddings.parquet"))).sort("compound_id")
tasks = checkpoint["tasks"]
state = {
    k.removeprefix("head."): v for k, v in checkpoint["state_dict"].items() if k.startswith("head.")
}
head = torch.nn.Sequential(
    torch.nn.Linear(600, 600),
    torch.nn.ReLU(),
    torch.nn.Dropout(0.0),
    torch.nn.Linear(600, state["3.weight"].shape[0]),
)
head.load_state_dict(state, strict=True)
head.eval()
torch.set_num_threads(4)
with torch.inference_mode():
    pred = head(
        torch.as_tensor(
            emb.select(sorted(c for c in emb.columns if c.startswith("emb_"))).to_numpy(),
            dtype=torch.float32,
        )
    ).numpy()
selected = {}
for task in tasks:
    if task["source"] in cfg.get("sources", ["pubchem1851", "tox21_p450_glo"]):
        # Ranking/continuous share one potency output; prefer physical scale when available.
        if task["column"] not in selected or task["objective"] == "continuous":
            selected[task["column"]] = task
columns, metadata = {}, []
for index, task in sorted(selected.items()):
    name = f"score_{index:04d}"
    values = pred[:, index]
    if task["objective"] == "continuous":
        values = values * task["std"] + task["mean"]
    assert np.isfinite(values).all()
    columns[name] = values
    metadata.append(dict(feature=name, **task))
assert columns
pl.DataFrame(columns).insert_column(0, emb["compound_id"]).write_parquet(
    out / "external-scores.parquet"
)
(out / "summary.json").write_text(
    json.dumps(
        dict(
            config=cfg,
            features=metadata,
            rows=emb.height,
            classification_scale="logit",
            selection="source identity only; no challenge labels",
        ),
        indent=2,
    )
)
