"""Standardize shared molecules without discarding assay or submission identities."""

import hashlib
import json
import math
import os
from collections import Counter
from importlib.metadata import version
from pathlib import Path

import polars as pl
from chembl_structure_pipeline import get_parent_mol, standardize_mol
from rdkit import Chem, rdBase

from .storage import write_json


def key(value):
    return hashlib.sha256(value.encode()).hexdigest()


def standardize(smiles):
    result = {
        "raw_smiles": smiles,
        "canonical_input": None,
        "standardized_smiles": None,
        "compound_id": None,
        "status": "invalid_input",
        "excluded": False,
        "error": None,
    }
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None or mol.GetNumAtoms() == 0:
            return result
        result["canonical_input"] = Chem.MolToSmiles(mol, isomericSmiles=True)
        result["status"] = "failed_standardize"
        standardized = standardize_mol(mol)
        result["status"] = "failed_parent"
        parent, excluded = get_parent_mol(standardized)
        result["excluded"] = bool(excluded)
        if excluded:
            result["status"] = "excluded"
            return result
        if parent is None or parent.GetNumAtoms() == 0:
            return result
        Chem.SanitizeMol(parent)
        canonical = Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True)
        result.update(
            standardized_smiles=canonical,
            compound_id=key(canonical),
            status="unchanged" if canonical == result["canonical_input"] else "changed",
        )
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
    return result


def clean(value):
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def prepare():
    inputs = Path(os.environ["CYP_INPUT_DIR"])
    output = Path(os.environ["CYP_OUTPUT_DIR"])
    config = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())
    if (
        not config["preserve_stereochemistry"]
        or config["tautomer_canonicalization"]
        or config["excluded_policy"] != "quarantine"
    ):
        raise ValueError("Unsupported standardization policy")
    files = sorted(inputs.glob("*.csv"))
    if len(files) != 5:
        raise ValueError("Expected the five official input CSVs")
    tables = {}
    for path in files:
        table = pl.read_csv(path, infer_schema_length=None, null_values=["", "NaN", "nan"])
        if not {"Molecule_Name", "SMILES"} <= set(table.columns):
            raise ValueError(f"Missing molecule columns: {path.name}")
        if table["SMILES"].null_count() or table["Molecule_Name"].null_count():
            raise ValueError("Null input molecule identity")
        tables[path.name.split("-", 1)[1]] = table
    smiles = sorted({s for t in tables.values() for s in t["SMILES"].to_list()})
    standardized = {}
    for i, smi in enumerate(smiles):
        standardized[smi] = standardize(smi)
        if (i + 1) % 500 == 0:
            print(f"Standardized {i + 1}/{len(smiles)}", flush=True)
    structures = pl.DataFrame(
        list(standardized.values()), infer_schema_length=None, schema_overrides={"error": pl.String}
    )
    structures.write_parquet(output / "standardization.parquet")
    compounds = (
        structures.filter(pl.col("compound_id").is_not_null())
        .select("compound_id", "standardized_smiles")
        .unique()
        .sort("compound_id")
    )
    compounds.write_parquet(output / "compounds.parquet")
    structures.filter(pl.col("compound_id").is_null()).write_parquet(output / "quarantine.parquet")
    records, origins, mappings = {}, [], []
    file_counts = {}
    with (output / "source_records.jsonl").open("x") as raw:
        for filename, table in tables.items():
            role = (
                "test"
                if "TEST-BLINDED" in filename
                else "aux"
                if "single-concentration" in filename
                else "train"
            )
            seen = set()
            for row_number, row in enumerate(table.iter_rows(named=True), start=1):
                row = {k: clean(v) for k, v in row.items()}
                payload = json.dumps(row, sort_keys=True, ensure_ascii=False, allow_nan=False)
                record_id = key(filename + "\n" + payload)
                origins.append(
                    {"source_file": filename, "source_row": row_number, "record_id": record_id}
                )
                if record_id in seen:
                    continue
                seen.add(record_id)
                std = standardized[row["SMILES"]]
                record = {
                    "record_id": record_id,
                    "source_file": filename,
                    "role": role,
                    "molecule_name": row["Molecule_Name"],
                    "raw_smiles": row["SMILES"],
                    "compound_id": std["compound_id"],
                    "payload": row,
                }
                records[record_id] = record
                raw.write(json.dumps(record, ensure_ascii=False, allow_nan=False) + "\n")
                mappings.append({k: record[k] for k in record if k != "payload"})
            file_counts[filename] = {
                "original_rows": table.height,
                "unique_records": len(seen),
                "exact_duplicates_removed": table.height - len(seen),
            }
    aliases = pl.DataFrame(mappings, infer_schema_length=None)
    aliases.write_parquet(output / "aliases.parquet")
    pl.DataFrame(origins).write_parquet(output / "record_origins.parquet")
    aliases.filter(pl.col("role") == "test").write_parquet(output / "test_mapping.parquet")
    # Measurement-level deduplication across overlapping DRC files, retaining provenance.
    measurements, measurement_sources = {}, []
    for record in records.values():
        if record["role"] != "train" or record["compound_id"] is None:
            continue
        row = record["payload"]
        for field, value in row.items():
            if (
                not field.startswith("CYP")
                or field.endswith(("_conf_high", "_conf_low", "_std"))
                or value is None
            ):
                continue
            measurement = {
                "compound_id": record["compound_id"],
                "endpoint": field,
                "value": float(value),
                "conf_low": row.get(field + "_conf_low"),
                "conf_high": row.get(field + "_conf_high"),
                "std": row.get(field + "_std"),
            }
            identifier = key(json.dumps(measurement, sort_keys=True))
            measurements[identifier] = {"measurement_id": identifier, **measurement}
            measurement_sources.append(
                {"measurement_id": identifier, "record_id": record["record_id"]}
            )
    values = pl.DataFrame(list(measurements.values()), infer_schema_length=None)
    values.write_parquet(output / "measurements.parquet")
    pl.DataFrame(measurement_sources).unique().write_parquet(output / "measurement_sources.parquet")
    conflicts = (
        values.group_by("compound_id", "endpoint")
        .agg(
            pl.len().alias("distinct_records"), pl.col("value").n_unique().alias("distinct_values")
        )
        .filter(pl.col("distinct_records") > 1)
        .sort("compound_id", "endpoint")
    )
    conflicts.write_parquet(output / "measurement_conflicts.parquet")
    duplicates = (
        aliases.filter(pl.col("compound_id").is_not_null())
        .group_by("compound_id")
        .agg(
            pl.col("molecule_name").unique().sort().alias("names"),
            pl.col("raw_smiles").n_unique().alias("raw_smiles_count"),
        )
        .filter((pl.col("names").list.len() > 1) | (pl.col("raw_smiles_count") > 1))
    )
    duplicates.write_parquet(output / "duplicate_compounds.parquet")
    overlap = (
        aliases.filter(pl.col("compound_id").is_not_null())
        .group_by("compound_id")
        .agg(pl.col("role").unique().sort().alias("roles"))
        .filter(pl.col("roles").list.contains("test") & (pl.col("roles").list.len() > 1))
    )
    overlap.write_parquet(output / "test_overlap.parquet")
    summary = {
        "versions": {
            "rdkit": rdBase.rdkitVersion,
            "chembl_structure_pipeline": version("chembl_structure_pipeline"),
        },
        "policy": config,
        "files": file_counts,
        "raw_unique_smiles": len(smiles),
        "unique_compounds": compounds.height,
        "status": dict(Counter(structures["status"])),
        "duplicate_compound_groups": duplicates.height,
        "deduplicated_measurements": values.height,
        "measurement_sources": len(measurement_sources),
        "conflicting_measurement_groups": conflicts.height,
        "test_overlap_compounds": overlap.height,
    }
    write_json(output / "summary.json", summary)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    prepare()
