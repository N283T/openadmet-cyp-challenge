"""Strict keyed-vector contract shared by extraction and DB registration."""

import numpy as np


def validate_embeddings(frame, expected_ids, columns):
    if frame.columns != ["compound_id", *columns]:
        raise ValueError("Embedding columns/order differ from the declared schema")
    ids = frame["compound_id"]
    if ids.null_count() or ids.n_unique() != frame.height:
        raise ValueError("Embedding IDs are missing or duplicated")
    if set(ids) != set(expected_ids) or len(expected_ids) != frame.height:
        raise ValueError("Embedding IDs do not cover the declared compound universe")
    values = frame.select(columns).to_numpy()
    if not np.isfinite(values).all():
        raise ValueError("Embedding contains missing or nonfinite values")
