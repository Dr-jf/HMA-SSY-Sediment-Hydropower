#!/usr/bin/env python3
"""Train the 17-predictor square-root SSY model. Author: Rui Deng.

The default training table has 1,963 records (1,911 + 52).
SSY is annual specific sediment yield in tonnes per square kilometre per year.
DeltaG and DeltaSC also accept the corresponding delta-symbol column names.
Requires NumPy, pandas, scikit-learn and XGBoost.
"""

import argparse
import csv
from dataclasses import dataclass
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import platform

import numpy as np
import pandas as pd


FEATURES = (
    "PGA", "S", "V", "P", "T", "L", "Ksn", "A", "G", "F", "R",
    "GMC", "Q", "SC", "SM", "DeltaG", "DeltaSC",
)
MONOTONIC_DIRECTIONS = (1, 1, -1, 1, 1, -1, 1, -1, 1, 1, 1, -1, 1, 1, 1, -1, -1)
MEMBER_TREES = (788, 736, 465, 316, 330)
MEMBER_FRACTION = 0.85
TEST_FRACTION = 0.30


@dataclass(frozen=True)
class SSYConfig:
    """Run settings for the five-member ensemble."""

    expected_samples: int = 1963
    split_seed: int = 42
    model_seed: int = 20261008
    member_trees: tuple[int, ...] = MEMBER_TREES
    threads: int = 1

    def __post_init__(self):
        if isinstance(self.expected_samples, bool) or not isinstance(self.expected_samples, int) or self.expected_samples < 3:
            raise ValueError("expected_samples must be an integer of at least 3.")
        for name in ("split_seed", "model_seed"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or not 0 <= value < 2**32:
                raise ValueError(f"{name} must be an integer in [0, 2**32).")
        if isinstance(self.threads, bool) or not isinstance(self.threads, int) or self.threads < 1:
            raise ValueError("threads must be a positive integer.")
        if len(self.member_trees) != 5 or any(
            isinstance(trees, bool) or not isinstance(trees, int) or trees < 1
            for trees in self.member_trees
        ):
            raise ValueError("member_trees must contain five positive integers.")


def _base_parameters():
    return {
        "max_depth": 6,
        "learning_rate": 0.02262021088229677,
        "min_child_weight": 0.06737187721626481,
        "subsample": 1.0,
        "colsample_bytree": 0.9779304618759297,
        "reg_lambda": 0.01032645739707662,
        "reg_alpha": 1.0,
        "gamma": 0.0019485455718767468,
        "max_bin": 512,
        "grow_policy": "depthwise",
    }


def _features(frame):
    if not isinstance(frame, pd.DataFrame):
        raise ValueError("Input must be a pandas DataFrame.")
    if frame.columns.duplicated().any():
        raise ValueError("Duplicate input column names are not allowed.")
    renamed = frame.copy(deep=False)
    aliases = {"\u0394G": "DeltaG", "\u0394SC": "DeltaSC"}
    for alias, public in aliases.items():
        if alias in frame.columns and public in frame.columns:
            raise ValueError(f"Input is ambiguous: provide either {public} or its delta-symbol alias.")
    renamed = renamed.rename(columns=aliases)
    missing = [name for name in FEATURES if name not in renamed.columns]
    if missing:
        raise ValueError("Missing required predictor columns: " + ", ".join(missing))
    try:
        values = renamed.loc[:, list(FEATURES)].apply(pd.to_numeric, errors="raise").to_numpy(
            dtype=np.float64, na_value=np.nan
        )
    except (ValueError, TypeError) as error:
        raise ValueError("Predictors must be numeric or missing (NaN).") from error
    if np.isinf(values).any():
        raise ValueError("Predictors cannot contain positive or negative infinity.")
    return pd.DataFrame(values, columns=FEATURES)


def _target(frame):
    if "SSY" not in frame.columns:
        raise ValueError("Missing required target column: SSY")
    try:
        values = pd.to_numeric(frame["SSY"], errors="raise").to_numpy(dtype=np.float64, na_value=np.nan)
    except (ValueError, TypeError) as error:
        raise ValueError("SSY must be numeric, finite and nonnegative.") from error
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("SSY must be numeric, finite and nonnegative.")
    return values


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def train_ssy(frame, output_dir, *, config=SSYConfig()):
    """Fit the training rows and save five model JSON files and metadata."""
    if not isinstance(config, SSYConfig):
        raise ValueError("config must be an SSYConfig instance.")
    X = _features(frame)
    y = _target(frame)
    if len(frame) != config.expected_samples:
        raise ValueError(f"Expected {config.expected_samples} samples; received {len(frame)}.")
    output = Path(output_dir)
    if output.is_symlink() or (output.exists() and (not output.is_dir() or any(output.iterdir()))):
        raise FileExistsError("Output must be a new or empty directory; existing files are preserved.")

    from sklearn.model_selection import train_test_split
    import xgboost as xgb

    training, held_out = train_test_split(
        np.arange(len(frame)), test_size=TEST_FRACTION,
        random_state=config.split_seed, shuffle=True,
    )
    output.mkdir(parents=True, exist_ok=True)
    members = []
    for index, trees in enumerate(config.member_trees):
        subset_seed = None if index == 0 else config.model_seed + 5000 + index
        if index == 0:
            selected = training
        else:
            count = max(2, int(np.floor(len(training) * MEMBER_FRACTION)))
            selected = training[np.random.default_rng(subset_seed).choice(len(training), count, replace=False)]
        model = xgb.XGBRegressor(
            **_base_parameters(), objective="reg:squarederror", tree_method="hist",
            monotone_constraints=MONOTONIC_DIRECTIONS, n_estimators=trees,
            random_state=config.model_seed, n_jobs=config.threads, verbosity=0,
        )
        # Fit each member on its selected training rows.
        model.fit(X.iloc[selected], np.sqrt(y[selected]), verbose=False)
        filename = f"member_{index + 1:02d}.json"
        path = output / filename
        model.save_model(path)
        members.append({
            "file": filename, "sha256": _sha256(path), "trees": trees,
            "training_samples": len(selected), "subset_seed": subset_seed,
        })
    metadata = {
        "format_version": 1,
        "author": "Rui Deng",
        "model": "fixed_17_predictor_sqrt_xgboost_five_member_core",
        "feature_names": list(FEATURES),
        "monotonic_directions": list(MONOTONIC_DIRECTIONS),
        "target": "SSY", "target_transform": "square_root",
        "prediction_rule": "maximum(mean(native_member_predictions), 0)**2",
        "base_parameters": _base_parameters(),
        "objective": "reg:squarederror", "tree_method": "hist",
        "model_seed": config.model_seed, "threads": config.threads,
        "member_subsample_fraction": MEMBER_FRACTION,
        "sample_counts": {"total": len(frame), "training": len(training), "test": len(held_out)},
        "expected_samples": config.expected_samples,
        "split": {"method": "shuffled row split", "test_size": TEST_FRACTION, "random_state": config.split_seed},
        "members": members,
        "versions": {"python": platform.python_version(), **{
            package: version(package) for package in ("numpy", "pandas", "scikit-learn", "xgboost")
        }},
        "recipe_status": "Fixed settings from model development, without retuning for the expanded input.",
    }
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return metadata


def predict_ssy(model_dir, frame):
    """Return raw SSY predictions using this core's saved five-member model."""
    X = _features(frame)
    root = Path(model_dir)
    metadata = json.loads((root / "metadata.json").read_text(encoding="utf-8"))
    if (
        metadata.get("format_version") != 1
        or metadata.get("feature_names") != list(FEATURES)
        or metadata.get("target_transform") != "square_root"
        or len(metadata.get("members", [])) != 5
    ):
        raise ValueError("Model metadata does not describe this 17-predictor core.")
    paths = []
    for index, member in enumerate(metadata["members"], 1):
        filename = f"member_{index:02d}.json"
        if member.get("file") != filename:
            raise ValueError("Unexpected model member filename.")
        path = root / filename
        if not path.is_file() or _sha256(path) != member.get("sha256"):
            raise ValueError("Model member hash verification failed.")
        paths.append(path)
    if len(X) == 0:
        return np.empty(0, dtype=float)
    import xgboost as xgb

    matrix = xgb.DMatrix(X, feature_names=list(FEATURES))
    total = np.zeros(len(X), dtype=float)
    for path in paths:
        model = xgb.Booster(params={"nthread": 1})
        model.load_model(path)
        if model.feature_names != list(FEATURES):
            raise ValueError("Saved model feature order differs from the required schema.")
        total += model.predict(matrix).astype(float) / 5
    prediction = np.maximum(total, 0.0) ** 2
    if not np.isfinite(prediction).all():
        raise ValueError("Model produced nonfinite SSY predictions.")
    return prediction


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Train the fixed 17-predictor SSY ensemble.",
        epilog="Required CSV columns: " + ", ".join((*FEATURES, "SSY")) + ".",
    )
    parser.add_argument("--input", type=Path, required=True, help="Private UTF-8 CSV with one station-year record per row.")
    parser.add_argument("--output", type=Path, required=True, help="New or empty output directory for the fitted models.")
    parser.add_argument("--expected-samples", type=int, default=1963, help="Required row count (default: 1963).")
    parser.add_argument("--split-seed", type=int, default=42, help="Random seed for the fixed 70:30 row split.")
    parser.add_argument("--threads", type=int, default=1, help="CPU threads per sequentially fitted member.")
    args = parser.parse_args(argv)
    try:
        config = SSYConfig(expected_samples=args.expected_samples, split_seed=args.split_seed, threads=args.threads)
        with args.input.open(encoding="utf-8-sig", newline="") as stream:
            header = next(csv.reader(stream), [])
        if len(header) != len(set(header)):
            raise ValueError("Duplicate CSV column names are not allowed.")
        frame = pd.read_csv(args.input, float_precision="round_trip")
        metadata = train_ssy(frame, args.output, config=config)
    except (ValueError, OSError) as error:
        parser.error(str(error))
    counts = metadata["sample_counts"]
    print(f"Saved five model members; training rows: {counts['training']}; held-out rows: {counts['test']}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
