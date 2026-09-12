"""Independent keyed membership, feature routing and saved-model inference audit."""

import json
import os
from pathlib import Path

import numpy as np
import polars as pl
import torch
from chemprop import featurizers
from chemprop.data import BatchMolGraph
from rdkit import Chem

from cyp_challenge.log2fc_pretrain import Network
from cyp_challenge.storage import write_json

inp, out = Path(os.environ["CYP_INPUT_DIR"]), Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
features = pl.read_parquet(inp / cfg["features"]).sort("compound_id")
original = pl.read_parquet(inp / cfg["nested"])
cf = original.filter(pl.col("fold") == pl.col("outer_fold")).sort("compound_id")
cols = [f"aux_{i}" for i in range(4)]
train = features.filter(pl.col("role") == "cross_fitted_train")
assert train["compound_id"].equals(cf["compound_id"])
assert np.array_equal(train.select(cols).to_numpy(), cf.select(cols).to_numpy())
actual = pl.read_parquet(inp / cfg["targets"])
members = pl.read_parquet(inp / cfg["membership"])
assert set(members["compound_id"]) == set(actual["compound_id"])
assert members.height == members["compound_id"].n_unique() == 4376
state = torch.load(inp / cfg["model"], map_location="cpu", weights_only=True)
metadata = state["metadata"]
recipe = {**metadata["recipe"], "initialization": "scratch"}
model = Network({"hyper_parameters": state["architecture"]}, recipe, 4)
model.load_state_dict(state["state_dict"], strict=True)
assert all(torch.isfinite(t).all() for t in model.state_dict().values())
torch.set_num_threads(3)
os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
torch.use_deterministic_algorithms(True)
torch.backends.cuda.matmul.allow_tf32 = False
torch.backends.cudnn.allow_tf32 = False
model.cuda().eval()
compounds = pl.read_parquet(inp / cfg["compounds"])
test = features.filter(pl.col("role") == "full_fit_test").join(
    compounds.select("compound_id", "standardized_smiles"),
    on="compound_id",
    validate="1:1",
    maintain_order="left",
)
assert test.height == 750 and not (set(test["compound_id"]) & set(members["compound_id"]))
f = featurizers.SimpleMoleculeMolGraphFeaturizer()
means = np.array(metadata["target_means"], dtype=np.float32)
stds = np.array(metadata["target_stds"], dtype=np.float32)
predictions = []
with torch.inference_mode():
    for start in range(0, test.height, recipe["batch_size"]):
        batch = test.slice(start, recipe["batch_size"])
        g = BatchMolGraph([f(Chem.MolFromSmiles(s)) for s in batch["standardized_smiles"]])
        g.to("cuda")
        predictions.append(model(g).cpu().numpy() * stds + means)
pred = np.concatenate(predictions)
diff = float(np.max(np.abs(pred - test.select(cols).to_numpy())))
assert diff < 1e-6
report = {
    "train_features_bitwise_equal": True,
    "train_rows": train.height,
    "all_actual_logfc_compounds_in_fit": True,
    "fit_compounds": members.height,
    "test_rows": test.height,
    "reloaded_model_max_absolute_difference": diff,
    "epochs": metadata["epochs"],
    "labels": metadata["train_labels"],
}
write_json(out / "verification.json", report)
print(json.dumps(report), flush=True)
