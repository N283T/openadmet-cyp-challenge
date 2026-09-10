import numpy as np

from cyp_challenge.log2fc_pretrain import partition, target_scaler


def test_outer_holdout_excluded_from_training_and_selection():
    ids = [str(i) for i in range(25)]
    folds = np.arange(25) % 5
    train, validation, heldout = partition(ids, folds, 2, 42, 0.1)
    assert set(heldout) == set(np.flatnonzero(folds == 2))
    assert not (set(train) & set(validation))
    assert not ((set(train) | set(validation)) & set(heldout))
    assert set(train) | set(validation) | set(heldout) == set(range(25))
    assert np.array_equal(train, partition(ids, folds, 2, 42, 0.1)[0])


def test_scaler_ignores_validation_and_outer_values():
    targets = np.array([[1.0, 5.0], [3.0, 5.0], [1e9, -1e9]])
    means, stds = target_scaler(targets, np.array([0, 1]))
    assert np.allclose(means, [2.0, 5.0])
    assert np.allclose(stds, [1.0, 1.0])


def test_nested_pair_excluded_from_fit_selection_and_scaling():
    ids = [str(i) for i in range(50)]
    folds = np.arange(50) % 5
    train, validation, heldout = partition(ids, folds, 1, 42, 0.1, excluded_folds=[1, 4])
    assert set(folds[heldout]) == {1, 4}
    assert set(folds[train]) | set(folds[validation]) == {0, 2, 3}
    values = np.arange(50, dtype=float).reshape(-1, 1)
    original = target_scaler(values, train)
    values[heldout] = 1e12
    values[validation] = -1e12
    changed = target_scaler(values, train)
    for before, after in zip(original, changed, strict=True):
        np.testing.assert_array_equal(before, after)
