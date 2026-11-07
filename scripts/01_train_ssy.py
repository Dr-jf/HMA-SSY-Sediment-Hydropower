#!/usr/bin/env python3
"""Train the fixed 17-predictor SSY core on a user-supplied private CSV."""

import argparse
import csv
from pathlib import Path
import sys

import pandas as pd

# Support a source checkout as well as an installed package.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from hma_models.ssy import SSYConfig, train_ssy


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Private input CSV; see docs/ssy.md for its schema.")
    parser.add_argument("--output", type=Path, required=True, help="New or empty local output directory; do not publish its weights.")
    parser.add_argument("--expected-samples", type=int, default=1963, help="Required row count (default: user-reported 1963).")
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
