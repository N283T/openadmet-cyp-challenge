"""Check pooling orientation, cofactor selection, and immutable-input storage."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

from cyp_challenge.runs import snapshot_input

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/pool_cyp_boltz_trunk.py"
spec = importlib.util.spec_from_file_location("pool_cyp", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_pool_chain_mapping_and_direction(tmp_path):
    n = 47
    s = np.zeros((n, 384), dtype=np.float16)
    z = np.zeros((n, n, 128), dtype=np.float16)
    z[0, 45:] = 2
    z[45:, 0] = 3
    z[44, 45:] = 7
    z[1:44, 45:] = 5
    path = tmp_path / "record.npz"
    np.savez(
        path,
        s=s,
        z=z,
        asym_id=np.array([0, 0] + [1] * 43 + [2] * 2),
        residue_index=np.array([0, 1] + [0] * 45),
        recycling_steps=1,
        embedding_stage="post",
    )
    result = module.pool_one(
        (
            path,
            "molecule",
            "train",
            {"length": 2, "core_residues": [1], "extended_residues": [1, 2]},
            42,
        )
    )
    assert result["z_core_pl_mean_000"] == 2
    assert result["z_core_lp_mean_000"] == 3
    assert result["z_fe_pl_mean_000"] == 7
    assert result["z_hem_other_pl_mean_000"] == 5
    assert result["ligand_tokens"] == 2


def test_input_storage_does_not_modify_source_and_survives_unlink(tmp_path):
    source = tmp_path / "source"
    source.write_bytes(b"retained")
    mode = source.stat().st_mode
    linked, copied = tmp_path / "linked", tmp_path / "copied"
    snapshot_input(source, linked, "hardlink")
    snapshot_input(source, copied)
    assert linked.stat().st_ino == source.stat().st_ino
    assert copied.stat().st_ino != source.stat().st_ino
    assert source.stat().st_mode == mode
    source.unlink()
    assert linked.read_bytes() == copied.read_bytes() == b"retained"
    with pytest.raises(FileExistsError):
        snapshot_input(linked, copied, "hardlink")


def test_prediction_export_accepts_both_readout_precisions():
    import polars as pl

    recovery_path = SCRIPT.parent / "recover_cyp_boltz_predictions.py"
    recovery_spec = importlib.util.spec_from_file_location("recovery", recovery_path)
    recovery = importlib.util.module_from_spec(recovery_spec)
    recovery_spec.loader.exec_module(recovery)
    single = np.array([0.3], dtype=np.float32)
    double = np.array([0.2], dtype=np.float64)
    result = recovery.combine_predictions(
        [
            pl.DataFrame({"model": ["tabpfn"], "prediction": single}),
            pl.DataFrame({"model": ["lgbm"], "prediction": double}),
        ]
    )
    assert result["prediction"].dtype == pl.Float64
    np.testing.assert_array_equal(result["prediction"].to_numpy(), [float(single[0]), double[0]])
