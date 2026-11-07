"""Artificial-data tests for the confidential-input SSY interface."""

import importlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest


FEATURES = [
    "PGA", "S", "V", "P", "T", "L", "Ksn", "A", "G", "F", "R",
    "GMC", "Q", "SC", "SM", "DeltaG", "DeltaSC",
]
REPOSITORY = Path(__file__).resolve().parents[1]


def artificial_frame(size=20):
    """Generate values in memory; these are not study observations."""
    rng = np.random.default_rng(7)
    frame = pd.DataFrame(rng.uniform(0.1, 2.0, (size, 17)), columns=FEATURES)
    frame["SSY"] = (2.0 + frame["P"] + 0.5 * frame["Q"]) ** 2
    return frame


def ssy_module():
    return importlib.import_module("hma_models.ssy")


def small_config(size=20):
    return ssy_module().SSYConfig(expected_samples=size, member_trees=(2, 2, 2, 2, 2))


def test_default_count_is_enforced_before_any_output_is_created(tmp_path):
    ssy = ssy_module()
    output = tmp_path / "must_not_exist"
    with pytest.raises(ValueError, match="1963"):
        ssy.train_ssy(artificial_frame(), output)
    assert not output.exists()


def test_default_1963_rows_are_split_dynamically_and_only_weights_are_saved(tmp_path):
    ssy = ssy_module()
    frame = artificial_frame(1963)
    frame["unused_private_identifier"] = "DO_NOT_EXPORT_THIS_IDENTIFIER"
    output = tmp_path / "artificial_models"
    config = ssy.SSYConfig(member_trees=(1, 1, 1, 1, 1))
    metadata = ssy.train_ssy(frame, output, config=config)
    assert metadata["sample_counts"] == {"total": 1963, "training": 1374, "test": 589}
    assert [member["training_samples"] for member in metadata["members"]] == [1374, 1167, 1167, 1167, 1167]
    assert metadata["split"]["random_state"] == 42
    assert len(list(output.iterdir())) == 6
    assert sorted(p.name for p in output.iterdir()) == [
        "member_01.json", "member_02.json", "member_03.json", "member_04.json",
        "member_05.json", "metadata.json",
    ]
    assert json.loads((output / "metadata.json").read_text()) == metadata
    for path in output.iterdir():
        assert "DO_NOT_EXPORT_THIS_IDENTIFIER" not in path.read_text()
    assert "predictions" not in metadata and "record_ids" not in metadata


def test_held_out_target_values_never_change_the_fitted_trees(tmp_path):
    ssy = ssy_module()
    frame = artificial_frame()
    changed = frame.copy()
    # The six held-out indices of sklearn's seed-42 70:30 split for 20 rows.
    changed.loc[[0, 17, 15, 1, 8, 5], "SSY"] = 1_000_000.0
    first = tmp_path / "first"
    second = tmp_path / "second"
    ssy.train_ssy(frame, first, config=small_config())
    ssy.train_ssy(changed, second, config=small_config())
    for name in [f"member_{i:02d}.json" for i in range(1, 6)]:
        assert (first / name).read_bytes() == (second / name).read_bytes()


def test_prediction_uses_native_member_mean_then_nonnegative_square(tmp_path):
    ssy = ssy_module()
    import xgboost as xgb

    frame = artificial_frame()
    frame.loc[2, "V"] = np.nan
    output = tmp_path / "models"
    ssy.train_ssy(frame, output, config=small_config())
    prediction_frame = frame[FEATURES[::-1]].copy()
    predicted = ssy.predict_ssy(output, prediction_frame)
    matrix = xgb.DMatrix(frame[FEATURES], feature_names=FEATURES)
    margins = []
    for index in range(1, 6):
        model = xgb.Booster()
        model.load_model(output / f"member_{index:02d}.json")
        margins.append(model.predict(matrix).astype(float))
    expected = np.maximum(np.mean(margins, axis=0), 0.0) ** 2
    np.testing.assert_allclose(predicted, expected, rtol=1e-12, atol=1e-12)
    assert predicted.shape == (20,)
    assert np.isfinite(predicted).all() and (predicted >= 0).all()


def test_native_delta_aliases_and_fixed_feature_order_preserve_predictions(tmp_path):
    ssy = ssy_module()
    frame = artificial_frame()
    alias_frame = frame.rename(columns={"DeltaG": "\u0394G", "DeltaSC": "\u0394SC"})
    output = tmp_path / "models"
    ssy.train_ssy(alias_frame, output, config=small_config())
    np.testing.assert_array_equal(
        ssy.predict_ssy(output, frame), ssy.predict_ssy(output, alias_frame)
    )
    ambiguous = frame.copy()
    ambiguous["\u0394G"] = ambiguous["DeltaG"]
    with pytest.raises(ValueError, match="ambiguous"):
        ssy.predict_ssy(output, ambiguous)


@pytest.mark.parametrize(
    "column,value",
    [("P", np.inf), ("V", "not-a-number"), ("SSY", np.nan),
     ("SSY", np.inf), ("SSY", -1.0)],
)
def test_bad_values_are_rejected_before_output_creation(tmp_path, column, value):
    ssy = ssy_module()
    frame = artificial_frame()
    if isinstance(value, str):
        frame[column] = frame[column].astype(object)
    frame.loc[1, column] = value
    output = tmp_path / "models"
    with pytest.raises(ValueError):
        ssy.train_ssy(frame, output, config=small_config())
    assert not output.exists()


@pytest.mark.parametrize("column", ["SSY", "PGA", "DeltaSC"])
def test_missing_columns_are_rejected(tmp_path, column):
    ssy = ssy_module()
    with pytest.raises(ValueError, match="Missing"):
        ssy.train_ssy(artificial_frame().drop(columns=column), tmp_path / "models", config=small_config())


def test_existing_nonempty_output_is_preserved(tmp_path):
    ssy = ssy_module()
    sentinel = tmp_path / "keep.txt"
    sentinel.write_text("Preserve me")
    with pytest.raises(FileExistsError):
        ssy.train_ssy(artificial_frame(), tmp_path, config=small_config())
    assert sentinel.read_text() == "Preserve me"
    assert sorted(p.name for p in tmp_path.iterdir()) == ["keep.txt"]


def test_model_integrity_is_checked_before_prediction(tmp_path):
    ssy = ssy_module()
    output = tmp_path / "models"
    ssy.train_ssy(artificial_frame(), output, config=small_config())
    with (output / "member_01.json").open("a") as stream:
        stream.write(" ")
    with pytest.raises(ValueError, match="hash"):
        ssy.predict_ssy(output, artificial_frame())


def test_module_import_does_not_train_or_import_xgboost():
    result = subprocess.run(
        [sys.executable, "-c", "import sys; import hma_models.ssy; assert 'xgboost' not in sys.modules"],
        cwd=REPOSITORY, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr


def test_cli_help_and_explicit_count_override(tmp_path):
    script = REPOSITORY / "scripts" / "01_train_ssy.py"
    help_result = subprocess.run([sys.executable, str(script), "--help"], capture_output=True, text=True)
    assert help_result.returncode == 0, help_result.stderr
    assert "--expected-samples" in help_result.stdout
    source = tmp_path / "artificial.csv"
    artificial_frame().to_csv(source, index=False)
    output = tmp_path / "artificial_models"
    result = subprocess.run(
        [sys.executable, str(script), "--input", str(source), "--output", str(output),
         "--expected-samples", "20", "--threads", "1"],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert json.loads((output / "metadata.json").read_text())["sample_counts"] == {
        "total": 20, "training": 14, "test": 6,
    }


def test_cli_rejects_duplicate_csv_header_before_pandas_can_rename_it(tmp_path):
    frame = artificial_frame()
    frame["copy_of_P"] = frame["P"]
    source = tmp_path / "duplicate.csv"
    source.write_text(frame.to_csv(index=False, lineterminator="\n").replace(",copy_of_P\n", ",P\n"))
    output = tmp_path / "models"
    result = subprocess.run(
        [sys.executable, str(REPOSITORY / "scripts" / "01_train_ssy.py"),
         "--input", str(source), "--output", str(output), "--expected-samples", "20"],
        capture_output=True, text=True,
    )
    assert result.returncode != 0
    assert "Duplicate" in result.stderr
    assert not output.exists()
