"""Compute all installed Python RDKit and Mordred 2D descriptors with explicit missingness."""

import json
import math
import multiprocessing as mp
import os
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from importlib.metadata import version
from pathlib import Path

import polars as pl
from mordred import Calculator, descriptors, is_missing
from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors

from .storage import write_json

CALCULATOR = None
RDKIT = None


def initialize():
    global CALCULATOR, RDKIT
    CALCULATOR = Calculator(descriptors, ignore_3D=True)
    RDKIT = list(Descriptors._descList)


def finite(value):
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("nonfinite")
    return number


def calculate(row):
    compound_id, smiles = row
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Standardized SMILES no longer parses: {compound_id}")
    rdkit_values, mordred_values, errors = {}, {}, {}
    for name, function in RDKIT:
        try:
            rdkit_values[name] = finite(function(mol))
        except Exception as exc:
            rdkit_values[name] = None
            errors[f"rdkit:{name}"] = f"{type(exc).__name__}: {exc}"
    result = CALCULATOR(mol)
    for descriptor, value in zip(CALCULATOR.descriptors, result, strict=True):
        name = str(descriptor)
        try:
            if is_missing(value):
                raise ValueError(str(value))
            mordred_values[name] = finite(value)
        except Exception as exc:
            mordred_values[name] = None
            errors[f"mordred:{name}"] = f"{type(exc).__name__}: {exc}"
    return {
        "compound_id": compound_id,
        "rdkit": rdkit_values,
        "mordred": mordred_values,
        "errors": errors,
    }


def main():
    output = Path(os.environ["CYP_OUTPUT_DIR"])
    inputs = list(Path(os.environ["CYP_INPUT_DIR"]).glob("*compounds.parquet"))
    if len(inputs) != 1:
        raise ValueError("Provide one compounds.parquet snapshot")
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    compounds = pl.read_parquet(inputs[0]).sort("compound_id")
    if compounds["compound_id"].n_unique() != compounds.height:
        raise ValueError("Duplicate compound IDs in descriptor input")
    initialize()
    names = {
        "rdkit": [name for name, _ in RDKIT],
        "mordred": [str(d) for d in CALCULATOR.descriptors],
    }
    metadata = {
        "names": names,
        "versions": {"rdkit": rdBase.rdkitVersion, "mordredcommunity": version("mordredcommunity")},
        "ignore_3D": True,
        "workers": config["workers"],
        "batch_size": config["batch_size"],
    }
    write_json(output / "descriptor_catalog.json", metadata)
    batch = []
    counts = {family: Counter() for family in names}
    batches, completed = 0, 0

    def save(rows, index):
        for family, columns in names.items():
            frame = pl.DataFrame(
                [{"compound_id": r["compound_id"], **r[family]} for r in rows],
                schema={"compound_id": pl.String, **dict.fromkeys(columns, pl.Float64)},
            )
            frame.write_parquet(output / f"{family}-{index:05d}.parquet")
            for column in columns:
                counts[family][column] += frame[column].null_count()
        errors = pl.DataFrame(
            [
                {"compound_id": r["compound_id"], "errors": json.dumps(r["errors"], sort_keys=True)}
                for r in rows
            ],
            schema={"compound_id": pl.String, "errors": pl.String},
        )
        errors.write_parquet(output / f"errors-{index:05d}.parquet")
        print(f"Saved batch {index}: {len(rows)} molecules", flush=True)

    # Spawn avoids inheriting Polars' thread pool across fork.
    with ProcessPoolExecutor(
        max_workers=config["workers"], mp_context=mp.get_context("spawn"), initializer=initialize
    ) as executor:
        for result in executor.map(calculate, compounds.iter_rows(), chunksize=8):
            batch.append(result)
            completed += 1
            if len(batch) == config["batch_size"]:
                save(batch, batches)
                batches += 1
                batch = []
        if batch:
            save(batch, batches)
            batches += 1
    summary = {
        "compounds": completed,
        "batches": batches,
        "descriptor_counts": {f: len(v) for f, v in names.items()},
        "missing_counts": {f: dict(v) for f, v in counts.items()},
        "all_missing": {f: [n for n, v in c.items() if v == completed] for f, c in counts.items()},
        "imputation": False,
        "feature_selection": False,
    }
    write_json(output / "summary.json", summary)
    print(f"Completed {completed} compounds", flush=True)


if __name__ == "__main__":
    main()
