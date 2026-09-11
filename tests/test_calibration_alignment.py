import numpy as np
import polars as pl
import pytest

from cyp_challenge.calibration_alignment import aligned_matrix


def test_shuffled_predictions_follow_reference_ids():
    reference = pl.DataFrame({"compound_id": ["b", "a", "c"], "fold": [1, 0, 1]})
    predictions = pl.DataFrame(
        {"compound_id": ["c", "b", "a"], "fold": [1, 1, 0], "mean": [30.0, 20.0, 10.0]}
    )
    np.testing.assert_array_equal(
        aligned_matrix(reference, predictions, ["mean"])[:, 0], [20.0, 10.0, 30.0]
    )
    with pytest.raises(ValueError, match="Missing"):
        aligned_matrix(reference, predictions.head(2), ["mean"])
