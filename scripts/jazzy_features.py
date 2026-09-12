"""Recorded standard Jazzy ETDG features, with explicit failure and charge policy."""

import concurrent.futures
import importlib.metadata
import json
import os
import signal
import time
from pathlib import Path

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["POLARS_MAX_THREADS"] = "2"

import numpy as np
import polars as pl
from jazzy.api import molecular_vector_from_smiles
from rdkit import Chem

KEYS = ["sdc", "sdx", "sa", "dga", "dgp", "dgtot"]


def timeout_handler(signum, frame):
    raise TimeoutError("Jazzy seed exceeded time limit")


def calculate(task):
    cid, smiles, seeds, timeout = task
    mol = Chem.MolFromSmiles(smiles)
    charge = Chem.GetFormalCharge(mol) if mol is not None else None
    rows = []
    for seed in seeds:
        row = dict(compound_id=cid, seed=seed, formal_charge=charge)
        start = time.monotonic()
        try:
            if mol is None:
                raise ValueError("invalid_smiles")
            if charge != 0:
                # Public Jazzy API fixes Kallisto total charge at zero. Do not
                # silently treat an ion as neutral or change the stored SMILES.
                raise ValueError("unsupported_net_charge_for_standard_jazzy")
            if len(Chem.GetMolFrags(mol)) != 1:
                raise ValueError("disconnected_structure")
            signal.signal(signal.SIGALRM, timeout_handler)
            signal.alarm(timeout)
            vec = molecular_vector_from_smiles(
                smiles,
                minimisation_method=None,
                embedding_type="3D",
                embedding_seed=seed,
                embedding_max_iterations=1000,
            )
            if not isinstance(vec, dict) or not all(k in vec for k in KEYS):
                raise ValueError(f"invalid_jazzy_result: {vec}")
            if not all(np.isfinite(vec[k]) for k in KEYS):
                raise ValueError("nonfinite_jazzy_result")
            row.update({k: float(vec[k]) for k in KEYS})
            row.update(status="ok", error="")
        except Exception as exc:
            row.update({k: None for k in KEYS})
            row.update(status="failed", error=f"{type(exc).__name__}: {exc}"[:500])
        finally:
            signal.alarm(0)
        row["seconds"] = time.monotonic() - start
        rows.append(row)
    return rows


def main():
    inp = Path(os.environ["CYP_INPUT_DIR"])
    out = Path(os.environ["CYP_OUTPUT_DIR"])
    cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    compounds = pl.read_parquet(next(inp.glob("000-*"))).sort("compound_id")
    if cfg.get("limit"):
        compounds = compounds.head(cfg["limit"])
    tasks = [
        (r["compound_id"], r["standardized_smiles"], cfg["seeds"], cfg["timeout_seconds"])
        for r in compounds.iter_rows(named=True)
    ]
    rows = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=cfg["workers"]) as pool:
        futures = [pool.submit(calculate, t) for t in tasks]
        for i, future in enumerate(concurrent.futures.as_completed(futures), 1):
            rows.extend(future.result())
            if i % 250 == 0:
                print(f"completed={i}/{len(tasks)}", flush=True)
    raw = pl.DataFrame(rows, infer_schema_length=None).with_columns(
        [pl.col(k).cast(pl.Float64) for k in KEYS]
    ).sort("compound_id", "seed")
    raw.write_parquet(out / "jazzy-seeds.parquet")
    features = raw.group_by("compound_id").agg(
        *[pl.col(k).mean().alias(f"jazzy_{k}_mean") for k in KEYS],
        *[pl.col(k).std(ddof=0).alias(f"jazzy_{k}_std") for k in KEYS],
        (pl.col("status") == "ok").sum().alias("jazzy_success_count"),
    ).sort("compound_id")
    assert features.height == compounds.height
    features.write_parquet(out / "jazzy-features.parquet")
    report = dict(
        config=cfg,
        versions={k: importlib.metadata.version(k) for k in ["jazzy", "rdkit", "kallisto"]},
        method="Jazzy standard ETDG, no force-field minimization, no 2D fallback",
        charge_policy="Net-charged and disconnected structures retained with missing Jazzy values",
        total_compounds=features.height,
        success_counts=features.group_by("jazzy_success_count").len().to_dicts(),
        errors=raw.filter(pl.col("status") != "ok").group_by("error").len().to_dicts(),
    )
    (out / "feature-report.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report), flush=True)


if __name__ == "__main__":
    main()
