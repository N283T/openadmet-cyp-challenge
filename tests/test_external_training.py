import numpy as np

from cyp_challenge.external_training import grouped_membership, snapshot_input


def test_snapshot_basename_does_not_match_external_suffix(tmp_path):
    for name in ["003-targets.parquet", "005-external-targets.parquet"]:
        (tmp_path / name).touch()
    assert snapshot_input(tmp_path, "targets.parquet").name == "003-targets.parquet"
    assert (
        snapshot_input(tmp_path, "external-targets.parquet").name == "005-external-targets.parquet"
    )


def test_stereoisomer_and_auxiliary_connections_never_leak_into_fitting():
    ids = [f"id{i}" for i in range(12)]
    keys = ["held", "held", "aux_held", "test", "a", "a", "b", "c", "d", "e", "f", "g"]
    folds = [0, 1, 2, 3, 1, 2, 1, 2, 3, 4, 1, 2]
    train, val, held = grouped_membership(
        ids, keys, folds, 0, 42, 0.3, outer_connectivity=["aux_held"], test_connectivity=["test"]
    )
    assert held.tolist() == [0]
    assert not ({0, 1, 2, 3} & (set(train) | set(val)))
    assert (4 in train) == (5 in train)
    assert (4 in val) == (5 in val)
    assert set(np.array(keys)[train]).isdisjoint(np.array(keys)[val])
    repeated = grouped_membership(
        ids, keys, folds, 0, 42, 0.3, outer_connectivity=["aux_held"], test_connectivity=["test"]
    )
    assert all(np.array_equal(a, b) for a, b in zip((train, val, held), repeated, strict=True))
