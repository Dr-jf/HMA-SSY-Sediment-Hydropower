#!/usr/bin/env python3
"""Historical reservoir sediment model (1980–2019). Author: Rui Deng.

Run with --help. The dam CSV supplies identifiers, construction years, total and
dead capacity (m³), basin area (km²), and semicolon-separated direct upstream IDs.
Discharge is daily in m³/s and basin SSY is annual in t km⁻² yr⁻¹. The model
starts at the study boundary and retains the legacy water deduction by deposited
volume. Python 3.11 or 3.12, NumPy, and pandas are required.
"""

from __future__ import annotations

import argparse
import csv
from collections import deque
from dataclasses import dataclass
from numbers import Real
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd


SECONDS_PER_DAY = 86400.0
DENSITY_T_M3 = 1.2
WINDOW_DAYS = 30


@dataclass(frozen=True)
class SedimentDam:
    """Prepared dam parameters; upstream_ids lists immediate upstream dams."""

    dam_id: str
    built_year: int
    total_capacity_m3: float
    dead_capacity_m3: float
    basin_area_km2: float
    upstream_ids: tuple[str, ...] = ()


def _integer(value: object, name: str) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be an integer")
    if not np.isfinite(value) or float(value) != int(value):
        raise ValueError(f"{name} must be a finite integer")
    return int(value)


def _dam_network(dams: Sequence[SedimentDam]):
    if not dams or any(not isinstance(dam, SedimentDam) for dam in dams):
        raise ValueError("Supply a nonempty sequence of SedimentDam records")
    for dam in dams:
        if not isinstance(dam.dam_id, str) or not dam.dam_id.strip() or dam.dam_id != dam.dam_id.strip():
            raise ValueError("dam_id must be a nonempty string without surrounding whitespace")
        built = _integer(dam.built_year, "built_year")
        if not 1 <= built <= 9999:
            raise ValueError("built_year must be a positive calendar year")
        values = (dam.total_capacity_m3, dam.dead_capacity_m3, dam.basin_area_km2)
        if any(isinstance(v, (bool, np.bool_)) or not isinstance(v, Real) for v in values):
            raise ValueError("Capacity and area parameters must be numeric")
        if not np.isfinite(values).all():
            raise ValueError("Capacity and area parameters must be finite")
        if dam.total_capacity_m3 <= 0 or dam.basin_area_km2 <= 0:
            raise ValueError("Total capacity and basin area must be positive")
        if not 0 <= dam.dead_capacity_m3 <= dam.total_capacity_m3:
            raise ValueError("Dead capacity must be between zero and total capacity")
        if isinstance(dam.upstream_ids, str) or not isinstance(dam.upstream_ids, (tuple, list)):
            raise ValueError("upstream_ids must be a tuple or list of identifiers")
        if any(not isinstance(v, str) or not v.strip() for v in dam.upstream_ids):
            raise ValueError("Upstream identifiers must be nonempty strings")
        if len(set(dam.upstream_ids)) != len(dam.upstream_ids):
            raise ValueError("Duplicate direct upstream identifiers are not allowed")
    ids = [dam.dam_id for dam in dams]
    if len(set(ids)) != len(ids):
        raise ValueError("Duplicate dam_id values are not allowed")
    ordered = sorted(dams, key=lambda dam: dam.dam_id)
    index = {dam.dam_id: i for i, dam in enumerate(ordered)}
    upstream: list[list[int]] = []
    downstream: list[list[int]] = [[] for _ in ordered]
    for i, dam in enumerate(ordered):
        if any(value not in index for value in dam.upstream_ids):
            raise ValueError(f"Unknown upstream identifier for {dam.dam_id}")
        items = [index[value] for value in dam.upstream_ids]
        upstream.append(items)
        for child in items:
            downstream[child].append(i)
    if any(len(items) > 1 for items in downstream):
        raise ValueError("A dam cannot feed multiple downstream dams without a split model")
    indegree = np.asarray([len(items) for items in upstream], dtype=int)
    ready = deque(i for i, degree in enumerate(indegree) if degree == 0)
    order = []
    while ready:
        i = ready.popleft()
        order.append(i)
        for parent in downstream[i]:
            indegree[parent] -= 1
            if indegree[parent] == 0:
                ready.append(parent)
    if len(order) != len(ordered):
        raise ValueError("The direct-upstream network contains a cycle")
    return ordered, upstream, order


def _columns(frame: pd.DataFrame, required: set[str], label: str) -> pd.DataFrame:
    if not isinstance(frame, pd.DataFrame) or not required.issubset(frame.columns):
        raise ValueError(f"{label} must contain columns: {', '.join(sorted(required))}")
    if frame.columns.duplicated().any():
        raise ValueError(f"{label} has duplicate column names")
    out = frame.copy()
    if not out["dam_id"].map(lambda v: isinstance(v, str) and bool(v.strip()) and v == v.strip()).all():
        raise ValueError(f"{label} dam_id values must be nonempty strings")
    return out


def _finite_nonnegative(series: pd.Series, label: str) -> np.ndarray:
    try:
        values = pd.to_numeric(series, errors="raise").to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be numeric") from exc
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError(f"{label} must be finite and nonnegative")
    return values


def _prepared_arrays(discharge, annual_ssy, dates, years, ids):
    q = _columns(discharge, {"date", "dam_id", "flow_m3s"}, "discharge")
    try:
        q["date"] = pd.to_datetime(q["date"], errors="raise")
        if q["date"].dt.tz is not None or not q["date"].eq(q["date"].dt.normalize()).all():
            raise ValueError("Discharge dates must be timezone-free calendar days")
    except (TypeError, ValueError, AttributeError) as exc:
        raise ValueError("Discharge dates must be timezone-free calendar days") from exc
    q["flow_m3s"] = _finite_nonnegative(q["flow_m3s"], "flow_m3s")
    q_keys = pd.MultiIndex.from_frame(q[["date", "dam_id"]])
    expected_q = pd.MultiIndex.from_product([dates, ids], names=["date", "dam_id"])
    if not q_keys.is_unique or len(q_keys) != len(expected_q) or not q_keys.isin(expected_q).all():
        raise ValueError("Discharge must contain exactly one row per dam and study calendar day")
    q_matrix = q.pivot(index="date", columns="dam_id", values="flow_m3s").reindex(index=dates, columns=ids).to_numpy(float)
    with np.errstate(over="ignore", invalid="ignore"):
        valid_daily_volume = np.isfinite(q_matrix * SECONDS_PER_DAY).all()
    if not valid_daily_volume:
        raise ValueError("Discharge magnitude overflows daily water volume")

    annual = _columns(annual_ssy, {"year", "dam_id", "ssy_t_km2_yr"}, "annual_ssy")
    annual_years = _finite_nonnegative(annual["year"], "SSY year")
    if np.any(annual_years != np.floor(annual_years)):
        raise ValueError("SSY year must be an integer")
    annual["year"] = annual_years.astype(int)
    annual["ssy_t_km2_yr"] = _finite_nonnegative(annual["ssy_t_km2_yr"], "ssy_t_km2_yr")
    annual_keys = pd.MultiIndex.from_frame(annual[["year", "dam_id"]])
    expected_annual = pd.MultiIndex.from_product([years, ids], names=["year", "dam_id"])
    if not annual_keys.is_unique or len(annual_keys) != len(expected_annual) or not annual_keys.isin(expected_annual).all():
        raise ValueError("Annual SSY must contain exactly one row per dam and study year")
    ssy = annual.pivot(index="year", columns="dam_id", values="ssy_t_km2_yr").reindex(index=years, columns=ids).to_numpy(float)
    return q_matrix, ssy


def simulate_historical_sediment(
    dams: Sequence[SedimentDam],
    discharge: pd.DataFrame,
    annual_ssy: pd.DataFrame,
    *,
    start_year: int = 1980,
    end_year: int = 2019,
) -> pd.DataFrame:
    """Simulate a full-calendar-year subset of 1980--2019, without prehistory.

    Inputs are prepared long tables. No source observations or spatial matching
    are bundled. Pre-commissioning trapping/state fields are missing, although
    incoming water and sediment pass through the network. Every output is a
    model estimate and all cumulative deposition starts at the study boundary.
    """
    start_year, end_year = _integer(start_year, "start_year"), _integer(end_year, "end_year")
    if not 1980 <= start_year <= end_year <= 2019:
        raise ValueError("Choose complete study years within 1980--2019")
    dams, upstream, order = _dam_network(dams)
    ids = [dam.dam_id for dam in dams]
    dates = pd.date_range(f"{start_year}-01-01", f"{end_year}-12-31", freq="D")
    years = np.arange(start_year, end_year + 1)
    q, ssy = _prepared_arrays(discharge, annual_ssy, dates, years, ids)
    n_days, n_dams = q.shape
    areas = np.asarray([dam.basin_area_km2 for dam in dams], dtype=float)
    capacities = np.asarray([dam.total_capacity_m3 for dam in dams], dtype=float)
    with np.errstate(over="ignore", invalid="ignore"):
        gross = ssy * areas[None, :]
        valid_capacity = np.isfinite(capacities * DENSITY_T_M3).all()
    if not np.isfinite(gross).all() or not valid_capacity:
        raise ValueError("Parameter magnitude overflows annual mass or capacity")
    local_annual = gross.copy()
    for dam, children in enumerate(upstream):
        if children:
            local_annual[:, dam] -= gross[:, children].sum(axis=1)
    clipped = local_annual < 0
    local_annual = np.maximum(local_annual, 0)
    local = np.zeros_like(q)
    uniform = np.zeros_like(local_annual, dtype=bool)
    for yi, year in enumerate(years):
        mask = dates.year == year
        for dam in range(n_dams):
            if local_annual[yi, dam] <= 0:
                continue
            flow = q[mask, dam]
            peak = float(flow.max())
            if peak > 0:
                # Scaling before squaring preserves Q^2 weights without overflow.
                weights = (flow / peak) ** 2
                local[mask, dam] = local_annual[yi, dam] * (weights / weights.sum())
            else:
                local[mask, dam] = local_annual[yi, dam] / int(mask.sum())
                uniform[yi, dam] = True

    # First STUDY-year Q seeds every dam's activation window; no pre-study data.
    warmup_volume = q[dates.year == start_year].mean(axis=0) * SECONDS_PER_DAY
    storage = np.zeros(n_dams)
    suspended = np.zeros(n_dams)
    capacity = capacities.copy()
    cumulative = np.zeros(n_dams)
    active = np.zeros(n_dams, dtype=bool)
    storage_history = np.zeros((n_dams, WINDOW_DAYS))
    outflow_history = np.zeros((n_dams, WINDOW_DAYS))
    positions = np.zeros(n_dams, dtype=int)
    numeric_fields = (
        "local_sediment_t", "incoming_sediment_t", "outgoing_sediment_t",
        "trapped_t", "storage_loss_m3", "cumulative_trapped_t",
        "remaining_capacity_m3", "stored_water_m3", "suspended_t",
        "trapping_efficiency", "residence_time_years", "water_in_m3",
        "water_out_m3", "initial_water_added_m3", "water_sink_m3",
        "water_balance_residual_m3", "sediment_balance_residual_t",
    )
    values = {name: np.full(n_days * n_dams, np.nan) for name in numeric_fields}
    operating = np.zeros(n_days * n_dams, dtype=bool)
    clipped_out = np.zeros(n_days * n_dams, dtype=bool)
    uniform_out = np.zeros(n_days * n_dams, dtype=bool)

    for day, date in enumerate(dates):
        exported_today = np.zeros(n_dams)
        yi = date.year - start_year
        for dam in order:
            row = day * n_dams + dam
            incoming_water = q[day, dam] * SECONDS_PER_DAY
            incoming_mass = local[day, dam] + sum(exported_today[child] for child in upstream[dam])
            values["local_sediment_t"][row] = local[day, dam]
            values["incoming_sediment_t"][row] = incoming_mass
            values["water_in_m3"][row] = incoming_water
            values["initial_water_added_m3"][row] = 0
            clipped_out[row] = clipped[yi, dam]
            uniform_out[row] = uniform[yi, dam]
            if date.year < dams[dam].built_year:
                exported_today[dam] = incoming_mass
                values["outgoing_sediment_t"][row] = incoming_mass
                values["water_out_m3"][row] = incoming_water
                values["water_sink_m3"][row] = 0
                values["water_balance_residual_m3"][row] = 0
                values["sediment_balance_residual_t"][row] = 0
                continue

            opening_water = storage[dam]
            opening_sediment = suspended[dam]
            initial_water = 0.0
            if not active[dam]:
                active[dam] = True
                storage[dam] = float(dams[dam].dead_capacity_m3)
                initial_water = storage[dam]
                storage_history[dam, :] = storage[dam]
                outflow_history[dam, :] = warmup_volume[dam]
            mean_storage = storage_history[dam].mean()
            mean_outflow = outflow_history[dam].mean()
            tau = mean_storage / (mean_outflow * 365.25) if mean_outflow > 0 else (1e9 if mean_storage > 0 else 0.0)
            efficiency = float(np.clip(tau / (0.012 + 1.02 * tau), 0, 1)) if capacity[dam] > 0 else 0.0
            deposited = min(efficiency * incoming_mass, capacity[dam] * DENSITY_T_M3)
            volume = deposited / DENSITY_T_M3
            capacity[dam] = max(capacity[dam] - volume, 0.0)

            # Retained research convention: deposition volume is also a water sink.
            # This is NOT a correction of physical water displacement or validation.
            available = max(storage[dam] + incoming_water - volume, 0.0)
            released = max(available - capacity[dam], 0.0)
            new_storage = min(capacity[dam], max(available - released, 0.0))
            mixed = max(suspended[dam] + incoming_mass - deposited, 0.0)
            denominator = storage[dam] + incoming_water
            concentration = mixed / denominator if denominator > 0 else 0.0
            outgoing_mass = min(mixed, concentration * released)
            new_suspended = max(mixed - outgoing_mass, 0.0)
            cumulative[dam] += deposited
            operating[row] = True
            record = {
                "outgoing_sediment_t": outgoing_mass,
                "trapped_t": deposited,
                "storage_loss_m3": volume,
                "cumulative_trapped_t": cumulative[dam],
                "remaining_capacity_m3": capacity[dam],
                "stored_water_m3": new_storage,
                "suspended_t": new_suspended,
                "trapping_efficiency": efficiency,
                "residence_time_years": tau,
                "water_out_m3": released,
                "initial_water_added_m3": initial_water,
                "water_sink_m3": volume,
                "water_balance_residual_m3": opening_water + initial_water + incoming_water - volume - released - new_storage,
                "sediment_balance_residual_t": opening_sediment + incoming_mass - deposited - outgoing_mass - new_suspended,
            }
            for name, value in record.items():
                values[name][row] = value
            storage[dam], suspended[dam] = new_storage, new_suspended
            exported_today[dam] = outgoing_mass
            position = positions[dam]
            storage_history[dam, position] = new_storage
            outflow_history[dam, position] = released
            positions[dam] = (position + 1) % WINDOW_DAYS

    result = pd.DataFrame({
        "date": np.repeat(dates.to_numpy(), n_dams),
        "year": np.repeat(dates.year.to_numpy(), n_dams),
        "dam_id": np.tile(np.asarray(ids, dtype=object), n_days),
        "is_operating": operating,
        **values,
        "negative_increment_clipped": clipped_out,
        "uniform_mass_fallback": uniform_out,
    })
    return result


def summarize_annual(daily: pd.DataFrame) -> pd.DataFrame:
    """Sum new daily deposition; retain study-period cumulative mass and NA rows."""
    required = {"date", "dam_id", "trapped_t", "storage_loss_m3", "cumulative_trapped_t", "is_operating"}
    table = _columns(daily, required, "daily sediment output")
    table["date"] = pd.to_datetime(table["date"], errors="raise")
    if table.duplicated(["date", "dam_id"]).any():
        raise ValueError("Daily sediment output has duplicate dam/date rows")
    table["year"] = table["date"].dt.year
    table = table.sort_values(["date", "dam_id"])
    grouped = table.groupby(["year", "dam_id"], sort=True)
    amounts = grouped[["trapped_t", "storage_loss_m3"]].sum(min_count=1)
    amounts["cumulative_trapped_t"] = grouped["cumulative_trapped_t"].last()
    amounts["is_operating"] = grouped["is_operating"].any()
    return amounts.reset_index()


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
