"""Submission contracts verified against the official Space on 2026-09-07."""

import polars as pl

IDENTIFIERS = ["SMILES", "Molecule_Name"]
ENDPOINTS = {
    "regression": [
        f"{cyp}_pIC50_direct_inhibition" for cyp in ("CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4")
    ],
    "classification": ["CYP2D6_is_TDI", "CYP3A4_is_TDI"],
}


def read_table(path):
    if path.suffix == ".parquet":
        return pl.read_parquet(path)
    if path.suffix == ".csv":
        return pl.read_csv(path)
    raise ValueError("Expected .csv or .parquet")


def validate(path, track, reference):
    data, test = read_table(path), read_table(reference)
    expected = IDENTIFIERS + ENDPOINTS[track]
    if set(data.columns) != set(expected):
        raise ValueError(f"Expected columns: {expected}")
    if data.height != 750 or test.height != 750:
        raise ValueError("Expected 750 test compounds")
    for frame in (data, test):
        if any(frame[col].null_count() for col in IDENTIFIERS):
            raise ValueError("Null identifier")
        if frame["Molecule_Name"].n_unique() != frame.height:
            raise ValueError("Duplicate molecule identifier")
    if (
        not data.select(IDENTIFIERS)
        .sort("Molecule_Name")
        .equals(test.select(IDENTIFIERS).sort("Molecule_Name"))
    ):
        raise ValueError("Molecule identifiers/SMILES do not match the reference")
    for col in ENDPOINTS[track]:
        values = data[col]
        if values.null_count():
            raise ValueError(f"Null prediction: {col}")
        if track == "regression":
            if not values.dtype.is_numeric() or not values.is_finite().all():
                raise ValueError(f"Non-finite/non-numeric prediction: {col}")
            if values.n_unique() <= 1 or values.std() < 0.01:
                raise ValueError(f"Constant/near-constant prediction: {col}")
        else:
            if values.dtype != pl.Boolean:
                raise ValueError(f"Classification predictions must be boolean: {col}")
            if values.n_unique() <= 1:
                raise ValueError(f"All predictions in one class: {col}")
    return {"track": track, "rows": data.height, "columns": data.columns, "valid": True}
