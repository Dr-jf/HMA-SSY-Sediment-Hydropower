#!/usr/bin/env python
"""Run prepared historical sediment inputs. Author: Rui Deng."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import pandas as pd

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from hma_models.sediment import SedimentDam, simulate_historical_sediment, summarize_annual


def _read_csv(path: Path, *, dtype: dict[str, type]) -> pd.DataFrame:
    """Validate the original header before pandas can rename duplicate fields."""
    with path.open(encoding="utf-8-sig", newline="") as stream:
        header = next(csv.reader(stream), [])
        if not header or any(not name.strip() for name in header):
            raise ValueError("CSV column headers must be nonempty")
        if len(header) != len(set(header)):
            raise ValueError("Duplicate CSV column headers are not allowed")
        stream.seek(0)
        return pd.read_csv(stream, dtype=dtype, keep_default_na=False)


def main() -> None:
    parser = argparse.ArgumentParser(description="Historical sediment simulation, without prehistory (1980--2019 only).")
    parser.add_argument("--dams", type=Path, required=True)
    parser.add_argument("--discharge", type=Path, required=True)
    parser.add_argument("--ssy", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--start-year", type=int, default=1980)
    parser.add_argument("--end-year", type=int, default=2019)
    args = parser.parse_args()
    try:
        if args.output.exists() and (not args.output.is_dir() or any(args.output.iterdir())):
            raise ValueError("Output must be a new or empty directory; existing files are never overwritten")
        table = _read_csv(args.dams, dtype={"dam_id": str, "upstream_ids": str})
        required = {"dam_id", "built_year", "total_capacity_m3", "dead_capacity_m3", "basin_area_km2", "upstream_ids"}
        if not required.issubset(table.columns):
            raise ValueError("Dam CSV is missing required prepared parameter columns")
        dams = []
        for row in table.itertuples(index=False):
            upstream_text = str(row.upstream_ids).strip()
            upstream = tuple(part.strip() for part in upstream_text.split(";")) if upstream_text else ()
            try:
                dams.append(SedimentDam(
                    str(row.dam_id), float(row.built_year), float(row.total_capacity_m3),
                    float(row.dead_capacity_m3), float(row.basin_area_km2), upstream,
                ))
            except (TypeError, ValueError) as exc:
                raise ValueError("Dam parameters must be complete numeric values") from exc
        discharge = _read_csv(args.discharge, dtype={"dam_id": str})
        ssy = _read_csv(args.ssy, dtype={"dam_id": str})
        daily = simulate_historical_sediment(
            dams, discharge, ssy, start_year=args.start_year, end_year=args.end_year,
        )
        annual = summarize_annual(daily)
        args.output.mkdir(parents=True, exist_ok=True)
        daily.to_csv(args.output / "daily_sediment.csv", index=False)
        annual.to_csv(args.output / "annual_sediment.csv", index=False)
    except (OSError, ValueError, TypeError) as exc:
        parser.error(str(exc))
    print(f"Wrote daily_sediment.csv and annual_sediment.csv for {args.start_year}--{args.end_year}")


if __name__ == "__main__":
    main()
