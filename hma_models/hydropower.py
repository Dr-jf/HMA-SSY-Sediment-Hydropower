"""Shared historical hydropower and annual active-capacity-update equations.

Only artificial tests are distributed. Prepared inputs remain the user's data.
The sediment scenario starts with zero accumulated sediment at study entry.
"""

from __future__ import annotations

import argparse
import calendar
import csv
from dataclasses import dataclass
from datetime import date, timedelta
import hashlib
import json
import math
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd


SECONDS_PER_DAY = 86400.0
SEDIMENT_DENSITY_T_M3 = 1.2
METHOD_VERSION = "historical-annual-active-capacity-v1"
DAILY_COLUMNS = [
    "date", "dam_id", "season", "flow_m3s", "active_capacity_start_m3",
    "active_capacity_end_m3", "storage_start_m3", "storage_end_m3",
    "routed_m3s", "turbine_m3s", "bypass_m3s", "release_m3s",
    "retained_excess_m3s", "retained_extra_m3s", "hydraulic_spill_m3",
    "capacity_update_spill_m3", "spill_m3", "head_mean_m", "energy_mwh",
    "cumulative_sediment_m3", "cumulative_active_loss_m3",
    "cumulative_dead_loss_m3", "unallocated_sediment_m3",
    "water_balance_residual_m3",
]


def _finite(value: object, name: str, minimum: float | None = None) -> float:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a finite number")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite number") from exc
    if not math.isfinite(number) or (minimum is not None and number < minimum):
        raise ValueError(f"{name} is nonfinite or below its allowed minimum")
    return number


def _integer(value: object, name: str) -> int:
    number = _finite(value, name)
    if number != int(number):
        raise ValueError(f"{name} must be an integer")
    return int(number)


def _identifier(value: object) -> str:
    if pd.isna(value) or not str(value).strip():
        raise ValueError("dam_id must not be missing or empty")
    return str(value).strip()


@dataclass(frozen=True)
class HydropowerDam:
    dam_id: str
    built_year: int
    active_capacity_m3: float
    dead_capacity_m3: float
    reservoir_area_km2: float
    dam_height_m: float
    normal_level_m: float
    dead_level_m: float
    intake_level_m: float
    installed_capacity_mw: float
    max_discharge_m3s: float
    wet_months: tuple[int, ...]
    efficiency: float = 0.9

    def __post_init__(self) -> None:
        object.__setattr__(self, "dam_id", _identifier(self.dam_id))
        built = _integer(self.built_year, "built_year")
        if not 1 <= built <= 9999:
            raise ValueError("built_year must be a valid Gregorian year")
        object.__setattr__(self, "built_year", built)
        positive = ["active_capacity_m3", "reservoir_area_km2", "dam_height_m",
                    "installed_capacity_mw", "max_discharge_m3s", "efficiency"]
        for field in positive:
            number = _finite(getattr(self, field), field)
            if number <= 0:
                raise ValueError(f"{field} must be positive")
            object.__setattr__(self, field, number)
        for field in ["normal_level_m", "dead_level_m", "intake_level_m"]:
            object.__setattr__(self, field, _finite(getattr(self, field), field))
        object.__setattr__(self, "dead_capacity_m3", _finite(self.dead_capacity_m3, "dead_capacity_m3", 0))
        if self.efficiency > 1:
            raise ValueError("efficiency must not exceed one")
        bottom = self.normal_level_m - self.dam_height_m
        if not bottom <= self.dead_level_m < self.normal_level_m:
            raise ValueError("Require estimated bottom <= dead level < normal level")
        months = tuple(_integer(m, "wet month") for m in self.wet_months)
        if (not months or len(set(months)) != len(months)
                or any(m < 1 or m > 12 for m in months)):
            raise ValueError("wet_months must be unique calendar months")
        months = tuple(sorted(months))
        if months != tuple(range(months[0], months[-1] + 1)) or len(months) == 12 or (1 in months and 12 in months):
            raise ValueError("Require one continuous, non-cross-year wet period")
        object.__setattr__(self, "wet_months", months)
        if not math.isfinite(self.total_capacity_m3) or not math.isfinite(self.shape_m):
            raise ValueError("Reservoir geometry must be finite")

    @property
    def total_capacity_m3(self) -> float:
        return self.active_capacity_m3 + self.dead_capacity_m3

    @property
    def shape_m(self) -> float:
        raw = self.reservoir_area_km2 * 1e6 * self.dam_height_m / self.total_capacity_m3
        return min(10.0, max(1.0, raw))

    @property
    def intake_used_m(self) -> float:
        return min(self.intake_level_m, self.dead_level_m)


def _mean_fraction_power(x0: float, x1: float, exponent: float) -> float:
    low, high = sorted((x0, x1))
    if high == 0:
        return 0.0
    relative_span = (high - low) / high
    if relative_span < 1e-12:
        return ((low + high) / 2) ** exponent
    if relative_span == 1:
        return high ** exponent / (exponent + 1)
    return high ** exponent * (-math.expm1((exponent + 1) * math.log1p(-relative_span))) / ((exponent + 1) * relative_span)


def daily_energy_mwh(
    dam: HydropowerDam, capacity_m3: float, storage_start_m3: float,
    storage_end_m3: float, discharge_m3s: float,
) -> tuple[float, float]:
    """Integrate instantaneous capped power; return energy and mean head.

    Discharge is the already scheduled constant flow. The instantaneous power
    cap does not trigger another release optimization.
    """
    capacity = _finite(capacity_m3, "capacity_m3")
    if capacity <= 0:
        raise ValueError("Positive capacity is required for storage-head integration")
    start = _finite(storage_start_m3, "storage_start_m3", 0)
    end = _finite(storage_end_m3, "storage_end_m3", 0)
    q = _finite(discharge_m3s, "discharge_m3s", 0)
    tolerance = max(1e-5, capacity * 1e-12)
    if max(start, end) > capacity + tolerance:
        raise ValueError("Storage exceeds active capacity")
    x0, x1 = min(1.0, start / capacity), min(1.0, end / capacity)
    exponent = 1.0 / dam.shape_m
    base = dam.dead_level_m - dam.intake_used_m
    height = dam.normal_level_m - dam.dead_level_m
    mean_head = base + height * _mean_fraction_power(x0, x1, exponent)
    if q == 0 or mean_head <= 0:
        return 0.0, max(0.0, mean_head)
    coefficient = 0.00981 * dam.efficiency * q
    low, high = sorted((x0, x1))
    low_head, high_head = base + height * low ** exponent, base + height * high ** exponent
    if coefficient * high_head <= dam.installed_capacity_mw:
        energy = 24 * coefficient * mean_head
    elif coefficient * low_head >= dam.installed_capacity_mw:
        energy = 24 * dam.installed_capacity_mw
    elif high - low < 1e-14 * max(1.0, high):
        energy = 24 * min(dam.installed_capacity_mw, coefficient * mean_head)
    else:
        cutoff = ((dam.installed_capacity_mw / coefficient - base) / height) ** dam.shape_m
        cutoff = min(high, max(low, cutoff))
        uncapped = (cutoff - low) * (base + height * _mean_fraction_power(low, cutoff, exponent))
        energy = 24 * (coefficient * uncapped + dam.installed_capacity_mw * (high - cutoff)) / (high - low)
    return min(24 * dam.installed_capacity_mw, max(0.0, energy)), mean_head


def allocate_cumulative_sediment(dam: HydropowerDam, volume_m3: float) -> dict[str, float]:
    """Allocate cumulative deposited volume once, never add cumulative values twice."""
    volume = _finite(volume_m3, "cumulative deposited volume", 0)
    bottom = dam.normal_level_m - dam.dam_height_m
    dead_fraction = (dam.dead_level_m - bottom) / (dam.normal_level_m - bottom)
    dead = min(dam.dead_capacity_m3, volume * dead_fraction)
    active = min(dam.active_capacity_m3, volume * (1 - dead_fraction))
    excess = max(0.0, volume - dead - active)
    transfer = min(dam.active_capacity_m3 - active, excess)
    active += transfer
    excess -= transfer
    transfer = min(dam.dead_capacity_m3 - dead, excess)
    dead += transfer
    excess -= transfer
    if abs(volume - dead - active - excess) > max(1e-5, volume * 1e-12):
        raise ArithmeticError("Sediment allocation does not conserve volume")
    return {"active_loss_m3": active, "dead_loss_m3": dead,
            "active_capacity_m3": max(0.0, dam.active_capacity_m3 - active),
            "unallocated_m3": max(0.0, excess)}


def _period(start_year: int, end_year: int) -> tuple[int, int]:
    start = _integer(start_year, "start_year")
    end = _integer(end_year, "end_year")
    if not 1980 <= start <= end <= 2019:
        raise ValueError("The historical study window must lie within 1980-2019")
    return start, end


def _dams_checked(dams: Sequence[HydropowerDam]) -> list[HydropowerDam]:
    selected = list(dams)
    if not selected or not all(isinstance(d, HydropowerDam) for d in selected):
        raise ValueError("Supply at least one HydropowerDam")
    ids = [d.dam_id for d in selected]
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate dam_id in parameter records")
    return selected


def _dated_frame(frame: pd.DataFrame, required: set[str]) -> pd.DataFrame:
    if not required.issubset(frame.columns):
        raise ValueError(f"Missing columns: {sorted(required - set(frame.columns))}")
    result = frame.copy()
    result["dam_id"] = result["dam_id"].map(_identifier)
    try:
        result["date"] = pd.to_datetime(result["date"], errors="raise")
    except (TypeError, ValueError) as exc:
        raise ValueError("Dates must be daily Gregorian calendar values") from exc
    dates = result["date"]
    if dates.isna().any() or dates.dt.tz is not None or not dates.equals(dates.dt.normalize()):
        raise ValueError("Dates must be nonmissing, timezone-free, midnight daily values")
    if ((dates.dt.year < 1980) | (dates.dt.year > 2019)).any():
        raise ValueError("Only historical dates within 1980-2019 are accepted")
    if result.duplicated(["dam_id", "date"]).any():
        raise ValueError("Duplicate dam_id/date record")
    return result


def _discharge_series(dams: list[HydropowerDam], frame: pd.DataFrame, start: int, end: int) -> dict[str, pd.Series]:
    data = _dated_frame(frame, {"dam_id", "date", "flow_m3s"})
    if set(data.dam_id) - {d.dam_id for d in dams}:
        raise ValueError("Discharge contains unknown dam_id records")
    result = {}
    for dam in dams:
        first = max(start, dam.built_year)
        if first > end:
            result[dam.dam_id] = pd.Series(dtype=float, index=pd.DatetimeIndex([]))
            continue
        expected = pd.date_range(f"{first}-01-01", f"{end}-12-31", freq="D")
        rows = data[(data.dam_id == dam.dam_id) & data.date.between(expected[0], expected[-1])].sort_values("date")
        if not pd.DatetimeIndex(rows.date).equals(expected):
            raise ValueError(f"Missing or discontinuous operational discharge for dam {dam.dam_id}")
        try:
            q = pd.to_numeric(rows.flow_m3s, errors="raise").to_numpy(dtype=float)
        except (TypeError, ValueError) as exc:
            raise ValueError("Operational discharge must be numeric") from exc
        if not np.isfinite(q).all() or (q < 0).any():
            raise ValueError("Operational discharge must be finite and nonnegative")
        result[dam.dam_id] = pd.Series(q, index=expected)
    return result


def _annual_increments(dams: list[HydropowerDam], frame: pd.DataFrame, start: int, end: int) -> dict[str, dict[int, float]]:
    required = {"dam_id", "year", "trapped_t"}
    if not required.issubset(frame.columns):
        raise ValueError("Annual sediment requires dam_id, year, trapped_t (annual incremental tonnes)")
    data = frame.copy()
    data["dam_id"] = data.dam_id.map(_identifier)
    data["year"] = data.year.map(lambda y: _integer(y, "sediment year"))
    if not data.year.between(1980, 2019).all():
        raise ValueError("Only annual increments within 1980-2019 are accepted")
    if data.duplicated(["dam_id", "year"]).any():
        raise ValueError("Duplicate dam_id/year annual sediment record")
    result = {}
    for dam in dams:
        years = list(range(max(start, dam.built_year), end + 1))
        rows = data[(data.dam_id == dam.dam_id) & data.year.isin(years)].sort_values("year")
        if rows.year.tolist() != years:
            raise ValueError(f"Missing operational annual sediment increment for dam {dam.dam_id}")
        masses = [_finite(v, "annual trapped_t", 0) for v in rows.trapped_t]
        result[dam.dam_id] = dict(zip(years, masses))
    return result


def _season_end(day: date, wet_months: tuple[int, ...]) -> date:
    wet = day.month in wet_months
    year, month = day.year, day.month
    for _ in range(12):
        year, month = (year + 1, 1) if month == 12 else (year, month + 1)
        if (month in wet_months) != wet:
            return date(year, month, 1)
    raise ValueError("Season classification must contain both wet and dry months")


def _wet_retention_target(capacity: float, storage: float, flows: np.ndarray, qmax: float) -> float:
    deficit = max(0.0, capacity - storage)
    bounded = np.minimum(flows, qmax)
    excess = float(np.maximum(flows - qmax, 0).sum()) * SECONDS_PER_DAY
    incoming = float(flows.sum()) * SECONDS_PER_DAY
    tolerance = max(1e-5, 1e-12 * max(1.0, capacity, incoming))
    if deficit <= excess + tolerance:
        return 0.0
    if incoming < deficit - tolerance:
        return float(bounded.max())
    remaining = min(max((deficit - excess) / SECONDS_PER_DAY, 0.0), float(bounded.sum()))
    ordered = np.sort(bounded)
    for index, available in enumerate(ordered):
        equal_share = remaining / (len(ordered) - index)
        if equal_share <= available:
            return max(0.0, float(equal_share))
        remaining -= float(available)
    return float(ordered[-1])


def _simulate_one(dam: HydropowerDam, flows: pd.Series, increments: dict[int, float] | None) -> list[dict]:
    storage = capacity = dam.active_capacity_m3
    cumulative_volume = 0.0
    geometry = allocate_cumulative_sediment(dam, 0)
    dates, discharge = flows.index, flows.to_numpy(dtype=float)
    rows = []
    wet_stop, extra_target = -1, 0.0
    for i, stamp in enumerate(dates):
        day, qin = stamp.date(), float(discharge[i])
        wet = day.month in dam.wet_months
        s0, c0 = storage, capacity
        release = retained_excess = retained_extra = 0.0
        if capacity == 0:
            route = min(qin, dam.max_discharge_m3s)
            hydraulic_spill = max(0.0, qin - route) * SECONDS_PER_DAY
            s_before_update = 0.0
            head = dam.normal_level_m - dam.intake_used_m
            energy = 24 * min(dam.installed_capacity_mw, 0.00981 * dam.efficiency * route * head)
        else:
            if wet:
                if i >= wet_stop:
                    remaining = (_season_end(day, dam.wet_months) - day).days
                    wet_stop = i + remaining
                    if wet_stop > len(dates):
                        raise ValueError("A complete wet period is required for refill planning")
                    extra_target = _wet_retention_target(capacity, storage, discharge[i:wet_stop], dam.max_discharge_m3s)
                base_flow = min(qin, dam.max_discharge_m3s)
                excess = max(0.0, qin - dam.max_discharge_m3s)
                deficit_rate = max(0.0, capacity - storage) / SECONDS_PER_DAY
                retained_excess = min(excess, deficit_rate)
                retained_extra = min(extra_target, base_flow, max(0.0, deficit_rate - retained_excess))
                route = base_flow - retained_extra
                hydraulic_spill = (excess - retained_excess) * SECONDS_PER_DAY
            else:
                remaining = (_season_end(day, dam.wet_months) - day).days
                release = min(storage / (remaining * SECONDS_PER_DAY), max(0.0, dam.max_discharge_m3s - qin))
                retained_excess = min(max(0.0, qin - dam.max_discharge_m3s), max(0.0, capacity - storage) / SECONDS_PER_DAY)
                route = min(qin, dam.max_discharge_m3s) + release
                hydraulic_spill = max(0.0, qin - dam.max_discharge_m3s - retained_excess) * SECONDS_PER_DAY
            s_before_update = storage + (qin - route) * SECONDS_PER_DAY - hydraulic_spill
            tolerance = max(1e-5, 1e-12 * max(1.0, capacity, qin * SECONDS_PER_DAY))
            if s_before_update < -tolerance or s_before_update > capacity + tolerance:
                raise ArithmeticError("Daily scheduling violates storage bounds")
            s_before_update = min(capacity, max(0.0, s_before_update))
            energy, head = daily_energy_mwh(dam, capacity, storage, s_before_update, route)
        storage = s_before_update
        capacity_spill = 0.0
        update_day = date(day.year, dam.wet_months[-1], calendar.monthrange(day.year, dam.wet_months[-1])[1])
        if increments is not None and day == update_day:
            # The full calendar-year increment is concentrated at this wet end.
            cumulative_volume = math.fsum([cumulative_volume, increments[day.year] / SEDIMENT_DENSITY_T_M3])
            geometry = allocate_cumulative_sediment(dam, cumulative_volume)
            capacity = geometry["active_capacity_m3"]
            if capacity > c0:
                raise ArithmeticError("Nonnegative sediment increased active capacity")
            storage = min(s_before_update, capacity)
            capacity_spill = s_before_update - storage
        spill = hydraulic_spill + capacity_spill
        residual = storage - s0 - qin * SECONDS_PER_DAY + route * SECONDS_PER_DAY + spill
        tolerance = max(1e-5, 1e-12 * max(1.0, dam.active_capacity_m3, qin * SECONDS_PER_DAY, spill))
        if (not all(math.isfinite(v) for v in (storage, energy, residual, route, spill))
                or abs(residual) > tolerance or not 0 <= storage <= capacity
                or route < 0 or route > dam.max_discharge_m3s + 1e-9 * max(1, dam.max_discharge_m3s)
                or spill < 0 or not 0 <= energy <= 24 * dam.installed_capacity_mw + 1e-8):
            raise ArithmeticError("Hydropower water, capacity, discharge or power invariant failed")
        rows.append(dict(
            date=stamp, dam_id=dam.dam_id, season="WET" if wet else "DRY", flow_m3s=qin,
            active_capacity_start_m3=c0, active_capacity_end_m3=capacity,
            storage_start_m3=s0, storage_end_m3=storage, routed_m3s=route,
            turbine_m3s=route if head > 0 else 0.0, bypass_m3s=route if head <= 0 else 0.0,
            release_m3s=release, retained_excess_m3s=retained_excess,
            retained_extra_m3s=retained_extra, hydraulic_spill_m3=hydraulic_spill,
            capacity_update_spill_m3=capacity_spill, spill_m3=spill,
            head_mean_m=head, energy_mwh=energy, cumulative_sediment_m3=cumulative_volume,
            cumulative_active_loss_m3=geometry["active_loss_m3"],
            cumulative_dead_loss_m3=geometry["dead_loss_m3"],
            unallocated_sediment_m3=geometry["unallocated_m3"], water_balance_residual_m3=residual,
        ))
    return rows


def _simulate(dams: Sequence[HydropowerDam], discharge: pd.DataFrame, annual_sediment: pd.DataFrame | None, start_year: int, end_year: int) -> pd.DataFrame:
    start, end = _period(start_year, end_year)
    selected = _dams_checked(dams)
    flows = _discharge_series(selected, discharge, start, end)
    increments = None if annual_sediment is None else _annual_increments(selected, annual_sediment, start, end)
    rows = []
    for dam in selected:
        rows.extend(_simulate_one(dam, flows[dam.dam_id], None if increments is None else increments[dam.dam_id]))
    result = pd.DataFrame(rows, columns=DAILY_COLUMNS)
    result["date"] = pd.to_datetime(result["date"])
    return result.sort_values(["dam_id", "date"]).reset_index(drop=True)


def simulate_baseline(dams: Sequence[HydropowerDam], discharge: pd.DataFrame, start_year: int = 1980, end_year: int = 2019) -> pd.DataFrame:
    """Simulate fixed initial capacity; initialize full storage at study entry."""
    return _simulate(dams, discharge, None, start_year, end_year)


def simulate_sediment(dams: Sequence[HydropowerDam], discharge: pd.DataFrame, annual_sediment: pd.DataFrame, start_year: int = 1980, end_year: int = 2019) -> pd.DataFrame:
    """Simulate annual wet-end capacity updates from incremental trapped_t."""
    return _simulate(dams, discharge, annual_sediment, start_year, end_year)


def calculate_net_loss(baseline: pd.DataFrame, sediment: pd.DataFrame) -> pd.DataFrame:
    """Strict keyed subtraction, preserving negative daily generation differences.

    This in-memory function checks keys and energy; callers must supply paired
    scenarios with the same forcing and parameters. The CLI additionally
    verifies their file fingerprints and study period.
    """
    required = {"date", "dam_id", "energy_mwh"}
    left, right = _dated_frame(baseline, required), _dated_frame(sediment, required)
    for frame in (left, right):
        values = pd.to_numeric(frame.energy_mwh, errors="raise").to_numpy(dtype=float)
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError("Scenario energy must be finite and nonnegative")
        frame["energy_mwh"] = values
    keys = ["dam_id", "date"]
    left_keys = set(map(tuple, left[keys].itertuples(index=False, name=None)))
    right_keys = set(map(tuple, right[keys].itertuples(index=False, name=None)))
    if left_keys != right_keys:
        raise ValueError("Baseline and sediment daily keys differ; no rows may be discarded")
    result = right.merge(left[keys + ["energy_mwh"]].rename(columns={"energy_mwh": "baseline_energy_mwh"}), on=keys, how="left", validate="one_to_one")
    result["net_loss_mwh"] = result.baseline_energy_mwh - result.energy_mwh
    return result.sort_values(keys).reset_index(drop=True)


def _read_prepared_csv(path: Path | str) -> pd.DataFrame:
    """Reject ambiguous raw headers and retain opaque identifier text.

    Inspect before pandas can rename duplicate/unnamed columns. Missing
    numeric cells remain empty strings and are validated only where used;
    this permits blank nonoperational annual sediment without changing IDs
    such as the literal string 'NA' into missing values.
    """
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        header = next(csv.reader(stream), None)
    if header is None:
        raise ValueError("CSV header is missing")
    normalized = [name.strip() for name in header]
    if not normalized or any(not name for name in normalized):
        raise ValueError("CSV header contains an empty column name")
    if len(normalized) != len(set(normalized)):
        raise ValueError("CSV header contains duplicate column names")
    return pd.read_csv(path, dtype={"dam_id": "string"}, keep_default_na=False, encoding="utf-8-sig")


def load_dams_csv(path: Path | str) -> list[HydropowerDam]:
    """Read the documented prepared parameter schema (wet_months uses ';')."""
    frame = _read_prepared_csv(path)
    names = set(HydropowerDam.__dataclass_fields__)
    required = names - {"efficiency"}
    if not required.issubset(frame.columns):
        raise ValueError(f"Missing dam parameter columns: {sorted(required - set(frame.columns))}")
    result = []
    for row in frame.to_dict(orient="records"):
        values = {key: row[key] for key in names if key in row}
        token = values["wet_months"]
        if pd.isna(token):
            raise ValueError("wet_months must be a semicolon-separated month list")
        values["wet_months"] = tuple(part.strip() for part in str(token).split(";"))
        result.append(HydropowerDam(**values))
    return _dams_checked(result)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _empty_output(path: Path) -> None:
    if path.exists() and (not path.is_dir() or any(path.iterdir())):
        raise FileExistsError("Output must be a new or empty directory; existing results are never overwritten")


def _annual_summary(daily: pd.DataFrame) -> pd.DataFrame:
    data = daily.copy()
    data["year"] = data.date.dt.year
    metrics = ["energy_mwh", "spill_m3", "capacity_update_spill_m3"]
    metrics += [key for key in ["baseline_energy_mwh", "net_loss_mwh"] if key in data]
    if data.empty:
        return pd.DataFrame(columns=["dam_id", "year", *metrics])
    annual = data.groupby(["dam_id", "year"], as_index=False)[metrics].sum()
    if "net_loss_mwh" in annual:
        annual["net_loss_fraction"] = annual.net_loss_mwh / annual.baseline_energy_mwh.replace(0, np.nan)
    return annual


def _write_outputs(output: Path, daily: pd.DataFrame, metadata: dict) -> None:
    _empty_output(output)
    output.mkdir(parents=True, exist_ok=True)
    daily_path = output / "daily_generation.csv"
    annual_path = output / "annual_generation.csv"
    daily.to_csv(daily_path, index=False, date_format="%Y-%m-%d")
    _annual_summary(daily).to_csv(annual_path, index=False)
    manifest = dict(metadata, status="COMPLETED", rows=len(daily),
                    daily_csv_sha256=_sha256(daily_path), annual_csv_sha256=_sha256(annual_path))
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False), encoding="utf-8")


def _parser(sediment: bool) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Historical hydropower with study-entry initialization (1980-2019 only).")
    parser.add_argument("--dams", required=True, type=Path)
    parser.add_argument("--discharge", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--start-year", type=int, default=1980)
    parser.add_argument("--end-year", type=int, default=2019)
    if sediment:
        parser.add_argument("--annual-sediment", required=True, type=Path)
        parser.add_argument("--baseline", required=True, type=Path, help="Baseline daily_generation.csv with adjacent manifest.json")
    return parser


def _cli(argv: Sequence[str] | None, sediment: bool) -> None:
    parser = _parser(sediment)
    args = parser.parse_args(argv)
    try:
        _empty_output(args.output)
        start, end = _period(args.start_year, args.end_year)
        hashes = {"dams_sha256": _sha256(args.dams), "discharge_sha256": _sha256(args.discharge)}
        dams = load_dams_csv(args.dams)
        q = _read_prepared_csv(args.discharge)
        metadata = dict(method_version=METHOD_VERSION, start_year=start, end_year=end,
                        initialization="full original active storage and zero sediment at study/commissioning entry",
                        inputs=hashes.copy(), implementation_sha256=_sha256(Path(__file__)))
        if sediment:
            baseline_manifest = json.loads(args.baseline.with_name("manifest.json").read_text(encoding="utf-8"))
            if (baseline_manifest.get("status") != "COMPLETED"
                    or baseline_manifest.get("model") != "historical_hydropower_baseline"
                    or baseline_manifest.get("method_version") != METHOD_VERSION
                    or baseline_manifest.get("implementation_sha256") != metadata["implementation_sha256"]
                    or baseline_manifest.get("start_year") != start
                    or baseline_manifest.get("end_year") != end
                    or baseline_manifest.get("inputs") != hashes
                    or baseline_manifest.get("daily_csv_sha256") != _sha256(args.baseline)):
                raise ValueError("Baseline provenance mismatch: require identical parameter/forcing SHA256, implementation, period and intact daily file")
            annual_hash = _sha256(args.annual_sediment)
            increments = _read_prepared_csv(args.annual_sediment)
            adjusted = simulate_sediment(dams, q, increments, start, end)
            baseline = _read_prepared_csv(args.baseline)
            daily = calculate_net_loss(baseline, adjusted)
            metadata["model"] = "historical_hydropower_annual_sediment"
            metadata["inputs"].update(annual_sediment_sha256=annual_hash, baseline_daily_sha256=_sha256(args.baseline))
            metadata["signed_net_loss"] = True
            if _sha256(args.annual_sediment) != annual_hash:
                raise ValueError("Annual sediment input changed during calculation")
        else:
            daily = simulate_baseline(dams, q, start, end)
            metadata["model"] = "historical_hydropower_baseline"
        if _sha256(args.dams) != hashes["dams_sha256"] or _sha256(args.discharge) != hashes["discharge_sha256"]:
            raise ValueError("Parameter or discharge input changed during calculation")
        _write_outputs(args.output, daily, metadata)
    except (ValueError, ArithmeticError, OSError, KeyError) as exc:
        parser.exit(2, f"Input or calculation error: {exc}\n")


def baseline_cli(argv: Sequence[str] | None = None) -> None:
    """Command-line adapter for the fixed-capacity historical baseline."""
    _cli(argv, sediment=False)


def sediment_cli(argv: Sequence[str] | None = None) -> None:
    """Command-line adapter for annual sediment updates and strict paired losses."""
    _cli(argv, sediment=True)
