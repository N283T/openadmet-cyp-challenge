"""Prepare challenge-disjoint ATTNSOM data and launch the pinned upstream trainer."""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import polars as pl
from rdkit import Chem

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
scratch = Path(os.environ["TMPDIR"]) / "attnsom-train"
scratch.mkdir(parents=True, exist_ok=False)
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def connectivity(molecule: Chem.Mol | None) -> str | None:
    if molecule is None:
        return None
    key = Chem.MolToInchiKey(molecule)
    return key.split("-", 1)[0] if key else None


compounds = pl.read_parquet(inp / cfg["compounds"])
forbidden = {
    connectivity(Chem.MolFromSmiles(smiles)) for smiles in compounds["standardized_smiles"]
}
if None in forbidden:
    raise ValueError("Could not create all challenge connectivity keys")

dataset_dir = scratch / "cyp_dataset"
dataset_dir.mkdir()
filter_rows = []
for isoform, filename in cfg["sdf_files"].items():
    supplier = Chem.SDMolSupplier(str(inp / filename), removeHs=False)
    writer = Chem.SDWriter(str(dataset_dir / f"{isoform}.sdf"))
    retained = 0
    excluded = 0
    for source_row, molecule in enumerate(supplier):
        if molecule is None:
            raise ValueError(f"Invalid {isoform} SDF record {source_row}")
        key = connectivity(molecule)
        if key in forbidden:
            excluded += 1
            continue
        primary_som = molecule.GetProp("PRIMARY_SOM")
        clean = Chem.Mol(molecule)
        for property_name in clean.GetPropNames(includePrivate=True, includeComputed=False):
            clean.ClearProp(property_name)
        clean.SetProp("_Name", f"{isoform}-{source_row}")
        clean.SetProp("PRIMARY_SOM", primary_som)
        writer.write(clean)
        retained += 1
    writer.close()
    filter_rows.append(
        {"isoform": isoform, "retained_instances": retained, "excluded_instances": excluded}
    )
pl.DataFrame(filter_rows).write_csv(out / "filtered-source-counts.csv")

source_dir = scratch / "source"
source_dir.mkdir()
for source_name, filename in cfg["source_files"].items():
    shutil.copy2(inp / filename, source_dir / source_name)

result_dir = out / "results"
checkpoint_dir = out / "checkpoints"
command = [
    sys.executable,
    str(source_dir / "main.py"),
    "--dataset_dir",
    str(dataset_dir),
    "--result_dir",
    str(result_dir),
    "--ckpt_dir",
    str(checkpoint_dir),
    "--save_ckpt",
    "--device",
    cfg["device"],
    "--n_splits",
    str(cfg["n_splits"]),
    "--inner_val_ratio",
    str(cfg["inner_val_ratio"]),
    "--max_epochs",
    str(cfg["max_epochs"]),
    "--batch_size",
    str(cfg["batch_size"]),
    "--seed",
    str(cfg["seed"]),
    "--hidden_size",
    str(cfg["hidden_size"]),
    "--num_layers",
    str(cfg["num_layers"]),
]
print("ATTNSOM command:", command, flush=True)
subprocess.run(command, cwd=source_dir, check=True)

summary = json.loads((result_dir / "summary.json").read_text())
summary["source_commit"] = cfg["source_commit"]
summary["challenge_connectivity_excluded"] = len(forbidden)
summary["filtered_source_counts"] = filter_rows
summary["training_command"] = command
(out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
