import polars as pl
import pytest

from cyp_challenge.embedding_validation import validate_embeddings


def test_embedding_identity_and_values():
    frame = pl.DataFrame({"compound_id": ["b", "a"], "fp0": [1.0, 2.0], "fp1": [3.0, 4.0]})
    validate_embeddings(frame, ["a", "b"], ["fp0", "fp1"])
    for bad in [
        frame.head(1),
        frame.with_columns(pl.lit("a").alias("compound_id")),
        frame.with_columns(pl.lit(float("nan")).alias("fp0")),
        frame.with_columns(pl.lit(None, dtype=pl.Float64).alias("fp0")),
        frame.select("compound_id", "fp1", "fp0"),
    ]:
        with pytest.raises(ValueError):
            validate_embeddings(bad, ["a", "b"], ["fp0", "fp1"])
