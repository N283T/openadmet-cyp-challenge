"""Audit the primary FDA/Instem CYP3A4 TDI SDF against challenge molecules."""

import json
import os
from collections import Counter
from pathlib import Path

import polars as pl
from rdkit import Chem

inp = Path(os.environ["CYP_INPUT_DIR"])
out = Path(os.environ["CYP_OUTPUT_DIR"])
cfg = json.loads(Path(os.environ["CYP_CONFIG_PATH"]).read_text())


def connectivity_key(molecule: Chem.Mol) -> str | None:
    try:
        molecule = Chem.RemoveHs(molecule)
        return Chem.MolToInchiKey(molecule).split("-")[0]
    except Exception:
        return None


def safe_prop(molecule: Chem.Mol, name: str) -> str:
    if not molecule.HasProp(name):
        return ""
    try:
        return molecule.GetProp(name)
    except UnicodeDecodeError:
        return "<non-utf8>"


compounds = pl.read_parquet(inp / cfg["compounds"])
labels = pl.read_parquet(inp / cfg["labels"])
test_ids = set(pl.read_parquet(inp / cfg["test_mapping"])["compound_id"].to_list())
challenge = (
    compounds.select("compound_id", "standardized_smiles")
    .with_columns(
        pl.Series(
            "connectivity",
            [
                connectivity_key(Chem.MolFromSmiles(value))
                for value in compounds["standardized_smiles"]
            ],
        ),
        pl.col("compound_id").is_in(test_ids).alias("is_test"),
    )
    .select("compound_id", "connectivity", "standardized_smiles", "is_test")
)
challenge_by_key = {row[1]: (row[0], row[3]) for row in challenge.iter_rows() if row[1] is not None}
a4_labels = dict(
    labels.filter(pl.col("endpoint") == "CYP3A4_is_TDI").select("compound_id", "value").iter_rows()
)

rows = []
invalid = 0
raw_values = Counter()
supplier = Chem.SDMolSupplier(str(inp / cfg["sdf"]), removeHs=False)
for source_row, molecule in enumerate(supplier):
    if molecule is None:
        invalid += 1
        continue
    raw = safe_prop(molecule, "CYP3A4 TDI").strip()
    raw_values[raw or "<missing>"] += 1
    if not raw:
        continue
    try:
        value = int(float(raw))
    except ValueError:
        value = None
    key = connectivity_key(molecule)
    match = challenge_by_key.get(key)
    compound_id = match[0] if match else None
    rows.append(
        {
            "source_row": source_row,
            "source_name": safe_prop(molecule, "Name/ID"),
            "pubchem_cid": safe_prop(molecule, "PubChem CID"),
            "raw_label": raw,
            "value": value,
            "connectivity": key,
            "canonical_smiles": Chem.MolToSmiles(
                Chem.RemoveHs(molecule), canonical=True, isomericSmiles=True
            ),
            "challenge_compound_id": compound_id,
            "challenge_is_test": bool(match[1]) if match else None,
            "challenge_value": int(a4_labels[compound_id]) if compound_id in a4_labels else None,
        }
    )

records = pl.DataFrame(rows)
records.write_parquet(out / "fda-cyp3a4-tdi-records.parquet")
valid = records.filter(pl.col("value").is_in([0, 1]))
deduplicated = valid.group_by("connectivity").agg(
    pl.len().alias("records"),
    pl.col("value").n_unique().alias("distinct_labels"),
    pl.col("value").first().alias("value"),
    pl.col("challenge_compound_id").drop_nulls().first().alias("challenge_compound_id"),
    pl.col("challenge_is_test").drop_nulls().first().alias("challenge_is_test"),
    pl.col("challenge_value").drop_nulls().first().alias("challenge_value"),
)
clean = deduplicated.filter(pl.col("distinct_labels") == 1)
overlap = clean.filter(pl.col("challenge_compound_id").is_not_null())
labeled_overlap = overlap.filter(pl.col("challenge_value").is_not_null())
conflicts = labeled_overlap.filter(pl.col("value") != pl.col("challenge_value"))

summary = {
    "sdf_molecules_with_tdi_property": records.height,
    "invalid_sdf_records": invalid,
    "raw_tdi_values": dict(sorted(raw_values.items())),
    "binary_tdi_records": valid.height,
    "unique_tdi_connectivities": deduplicated.height,
    "duplicate_connectivities": deduplicated.filter(pl.col("records") > 1).height,
    "conflicting_duplicate_connectivities": deduplicated.filter(
        pl.col("distinct_labels") > 1
    ).height,
    "clean_unique_tdi_connectivities": clean.height,
    "positive_fraction": float(clean["value"].mean()),
    "challenge_overlap": overlap.height,
    "challenge_train_overlap": overlap.filter(~pl.col("challenge_is_test")).height,
    "challenge_test_overlap": overlap.filter(pl.col("challenge_is_test")).height,
    "challenge_labeled_overlap": labeled_overlap.height,
    "challenge_label_conflicts": conflicts.height,
    "challenge_label_agreement": (
        float(1 - conflicts.height / labeled_overlap.height) if labeled_overlap.height else None
    ),
    "interpretation": (
        "Source labels are heterogeneous paper-defined CYP3A4 TDI outcomes and are not "
        "assumed interchangeable with challenge labels. Exact connectivity overlap is audited "
        "before any possible auxiliary-training experiment."
    ),
}
(out / "audit-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2), flush=True)
