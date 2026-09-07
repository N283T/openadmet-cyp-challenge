import polars as pl
import pytest

from cyp_challenge.storage import allocate, write_json
from cyp_challenge.tracks import ENDPOINTS, validate


def files(tmp_path, track):
    frame = pl.DataFrame({"SMILES": ["C"] * 750, "Molecule_Name": [f"m{i}" for i in range(750)]})
    reference = tmp_path / "test.csv"
    frame.write_csv(reference)
    for col in ENDPOINTS[track]:
        values = (
            [bool(i % 2) for i in range(750)] if track == "classification" else list(range(750))
        )
        frame = frame.with_columns(pl.Series(col, values))
    prediction = tmp_path / "predictions.parquet"
    frame.write_parquet(prediction)
    return frame, prediction, reference


@pytest.mark.parametrize("track", ENDPOINTS)
def test_valid_and_reordered_predictions(tmp_path, track):
    frame, prediction, reference = files(tmp_path, track)
    frame.reverse().write_parquet(prediction)
    assert validate(prediction, track, reference)["valid"]


def test_id_mismatch(tmp_path):
    frame, prediction, reference = files(tmp_path, "regression")
    frame.with_columns(pl.lit("wrong").alias("SMILES")).write_parquet(prediction)
    with pytest.raises(ValueError, match="do not match"):
        validate(prediction, "regression", reference)


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), None])
def test_invalid_regression_values(tmp_path, bad):
    frame, prediction, reference = files(tmp_path, "regression")
    frame.with_columns(pl.lit(bad).alias(ENDPOINTS["regression"][0])).write_parquet(prediction)
    with pytest.raises(ValueError):
        validate(prediction, "regression", reference)


def test_probability_is_not_classification_submission(tmp_path):
    frame, prediction, reference = files(tmp_path, "classification")
    frame.with_columns(pl.lit(0.5).alias(ENDPOINTS["classification"][0])).write_parquet(prediction)
    with pytest.raises(ValueError, match="boolean"):
        validate(prediction, "classification", reference)


def test_constant_predictions_rejected(tmp_path):
    frame, prediction, reference = files(tmp_path, "regression")
    frame.with_columns(pl.lit(5.0).alias(ENDPOINTS["regression"][0])).write_parquet(prediction)
    with pytest.raises(ValueError, match="constant"):
        validate(prediction, "regression", reference)


def test_no_record_overwrite(tmp_path):
    path = tmp_path / "record.json"
    write_json(path, {"original": True})
    with pytest.raises(FileExistsError):
        write_json(path, {"original": False})
    assert '"original": true' in path.read_text()


def test_unique_directories(monkeypatch, tmp_path):
    monkeypatch.setattr("cyp_challenge.storage.ROOT", tmp_path)
    first, a = allocate("runs", "regression")
    second, b = allocate("runs", "regression")
    assert first != second and a != b
    assert a.is_dir() and b.is_dir()
