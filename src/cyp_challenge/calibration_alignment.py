"""Order explicit prediction keys before any positional numeric operations."""

import numpy as np


def aligned_matrix(reference, predictions, columns):
    keys = ["compound_id", "fold"]
    joined = reference.select(keys).join(
        predictions.select(*keys, *columns),
        on=keys,
        how="left",
        validate="1:1",
        maintain_order="left",
    )
    if joined.select(keys).rows() != reference.select(keys).rows():
        raise ValueError("Prediction keys do not preserve reference order")
    values = joined.select(columns).to_numpy().astype(float)
    if not np.isfinite(values).all():
        raise ValueError("Missing or nonfinite prediction after keyed alignment")
    return values
