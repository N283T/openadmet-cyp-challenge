"""Replay keyed frozen and auxiliary features and resolve the official Emax delta."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
from rdkit import Chem

from cyp_challenge.chemeleon import Fingerprint
from cyp_challenge.log2fc_pretrain import Network

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def p(key):
    return inp / cfg[key]


torch.set_num_threads(4)
torch.backends.cuda.matmul.allow_tf32 = False
torch.backends.cudnn.allow_tf32 = False
local, remote = (
    pl.read_csv(p("local"), infer_schema_length=None),
    pl.read_csv(p("remote"), infer_schema_length=None),
)
removed = local.join(remote, on=local.columns, how="anti", nulls_equal=True)
added = remote.join(local, on=local.columns, how="anti", nulls_equal=True)
removed.write_csv(out / "removed-emax.csv")
added.write_csv(out / "added-emax.csv")
compounds = pl.read_parquet(p("compounds")).sort("compound_id")
saved_std = pl.read_parquet(p("standardization"))
changed = saved_std.filter(pl.col("status") == "changed")
rows = []
for row in changed.iter_rows(named=True):
    a, b = [Chem.MolFromSmiles(row[c]) for c in ["raw_smiles", "standardized_smiles"]]
    rows.append(
        dict(
            compound_id=row["compound_id"],
            raw_smiles=row["raw_smiles"],
            standardized_smiles=row["standardized_smiles"],
            raw_stereocenters=len(Chem.FindMolChiralCenters(a)),
            parent_stereocenters=len(Chem.FindMolChiralCenters(b)),
            raw_charge=Chem.GetFormalCharge(a),
            parent_charge=Chem.GetFormalCharge(b),
        )
    )
pl.DataFrame(rows).write_csv(out / "structural-changes.csv")
features = pl.read_parquet(p("embeddings"))
assert features["compound_id"].n_unique() == compounds.height
assert set(features["compound_id"]) == set(compounds["compound_id"])
sample = (
    pl.concat(
        [
            compounds.sample(n=128, seed=42),
            compounds.join(changed.select("compound_id"), on="compound_id"),
        ]
    )
    .unique()
    .sort("compound_id")
)
encoder = Fingerprint(p("checkpoint"), "cuda")
actual = encoder(sample["standardized_smiles"].to_list())
cols = sorted(c for c in features.columns if c != "compound_id")
expected = (
    sample.select("compound_id")
    .join(features, on="compound_id", validate="1:1", maintain_order="left")
    .select(cols)
    .to_numpy()
)
np.testing.assert_allclose(actual, expected, atol=2e-5, rtol=2e-5)
metadata = json.loads(p("metadata").read_text())
checkpoint = torch.load(p("checkpoint"), map_location="cpu", weights_only=True)
model = Network(checkpoint, metadata["config"], len(metadata["heads"])).cuda().eval()
best = torch.load(p("best"), map_location="cpu", weights_only=False)
model.load_state_dict(best["state_dict"], strict=True)
with torch.inference_mode():
    raw = model(encoder.graph(sample["standardized_smiles"].to_list())).cpu().numpy()
actual_aux = raw * np.asarray(metadata["target_stds"], dtype=np.float32) + np.asarray(
    metadata["target_means"], dtype=np.float32
)
preds = pl.concat(pl.read_parquet(inp / f) for f in cfg["predictions"])
expected_aux = (
    sample.select("compound_id")
    .join(preds, on="compound_id", validate="1:1", maintain_order="left")
    .select(metadata["heads"])
    .to_numpy()
)
np.testing.assert_allclose(actual_aux, expected_aux, atol=2e-5, rtol=2e-5)
summary = dict(
    emax_removed=removed.to_dicts(),
    emax_added=added.height,
    excluded_smiles=saved_std.filter(pl.col("status") == "excluded").to_dicts(),
    frozen_all_ids_complete=True,
    replay_sample=sample.height,
    frozen_max_abs=float(np.abs(actual - expected).max()),
    auxiliary_max_abs=float(np.abs(actual_aux - expected_aux).max()),
    auxiliary_heads=metadata["heads"],
    auxiliary_inverse_transform_pass=True,
)
(out / "summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
