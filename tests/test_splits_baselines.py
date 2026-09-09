import importlib.util
from pathlib import Path

import numpy as np
import pytest

from cyp_challenge.splits import grouped_folds

spec = importlib.util.spec_from_file_location(
    "split_baselines", Path(__file__).parents[1] / "scripts/split_baselines.py"
)
baseline = importlib.util.module_from_spec(spec)
spec.loader.exec_module(baseline)


def test_group_integrity_and_empty_scaffolds():
    groups = ["", "", "a", "a", "b", "c", "d", "e"]
    folds = grouped_folds(groups, 3)
    assert set(folds) == {0, 1, 2}
    for group in set(groups):
        assert len(set(folds[np.array(groups) == group])) == 1
    assert np.array_equal(folds, grouped_folds(groups, 3))
    with pytest.raises(ValueError):
        grouped_folds(["", ""], 2)


def test_st_rae_interval_and_baseline():
    y = np.array([0.0, 2.0, 4.0])
    lower, upper = y - 0.25, y + 0.25
    assert baseline.st_rae(y, y + 0.1, lower, upper) == 0
    assert baseline.st_rae(y, np.full(3, y.mean()), lower, upper) == 1
    p = np.array([1.0, 2.0, 2.0])
    assert baseline.st_rae(y, p, y, y) == pytest.approx(3 / 4)
    assert baseline.st_rae(y, p, np.full(3, -10.0), np.full(3, 10.0)) is None


def test_tanimoto_not_bit_agreement():
    query = np.array([[1, 1, 0], [0, 0, 0]])
    reference = np.array([[1, 0, 0], [0, 0, 1]])
    assert np.allclose(baseline.nearest(query, reference), [0.5, 0])


def test_mixed_track_metric_schema():
    rows = [{"ST-RAE": 0.8}] * 144 + [{"MCC": 0.2}]
    frame = baseline.metric_table(rows)
    assert frame["MCC"][-1] == 0.2
