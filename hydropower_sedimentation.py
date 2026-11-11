"""Annual wet-season capacity updates and signed historical generation loss."""

from __future__ import annotations

import calendar
from datetime import date
import hashlib
import json
import math
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from hydropower_baseline import (
    METHOD_VERSION,
    HydropowerDam,
    _check_cli_inputs,
    _dams_checked,
    _dated_frame,
    _finite,
    _identifier,
    _integer,
    _load_cli_inputs,
    _parser,
    _period,
    _read_prepared_csv,
    _sha256,
    _simulate,
    _write_outputs,
)

SEDIMENT_DENSITY_T_M3 = 1.2


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


def simulate_sediment(
    dams: Sequence[HydropowerDam], discharge: pd.DataFrame, annual_sediment: pd.DataFrame,
    start_year: int = 1980, end_year: int = 2019,
) -> pd.DataFrame:
    """Schedule annual trapped sediment after generation on the last wet day."""
    start, end = _period(start_year, end_year)
    selected = _dams_checked(dams)
    increments = _annual_increments(selected, annual_sediment, start, end)
    schedules = {}
    for dam in selected:
        updates = {}
        cumulative_volume = 0.0
        month = dam.wet_months[-1]
        for year, trapped_t in increments[dam.dam_id].items():
            cumulative_volume = math.fsum([cumulative_volume, trapped_t / SEDIMENT_DENSITY_T_M3])
            update_day = date(year, month, calendar.monthrange(year, month)[1])
            geometry = allocate_cumulative_sediment(dam, cumulative_volume)
            geometry["cumulative_sediment_m3"] = cumulative_volume
            updates[update_day] = geometry
        schedules[dam.dam_id] = updates
    return _simulate(selected, discharge, start, end, schedules)


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


def sediment_cli(argv: Sequence[str] | None = None) -> None:
    parser = _parser("Historical hydropower with annual sediment capacity updates (1980-2019).")
    parser.add_argument("--annual-sediment", required=True, type=Path)
    parser.add_argument("--baseline", required=True, type=Path,
                        help="Baseline daily_generation.csv with adjacent manifest.json")
    args = parser.parse_args(argv)
    try:
        dams, discharge, metadata = _load_cli_inputs(args)
        start, end = metadata["start_year"], metadata["end_year"]
        hashes = metadata["inputs"]
        baseline_hash = _sha256(args.baseline)
        manifest = json.loads(args.baseline.with_name("manifest.json").read_text(encoding="utf-8"))
        if (manifest.get("status") != "COMPLETED"
                or manifest.get("model") != "historical_hydropower_baseline"
                or manifest.get("method_version") != METHOD_VERSION
                or manifest.get("implementation_sha256") != metadata["implementation_sha256"]
                or manifest.get("start_year") != start
                or manifest.get("end_year") != end
                or manifest.get("inputs") != hashes
                or manifest.get("daily_csv_sha256") != baseline_hash):
            raise ValueError(
                "Baseline provenance mismatch: require identical parameter/forcing SHA256, "
                "implementation, period and intact daily file"
            )
        annual_hash = _sha256(args.annual_sediment)
        annual = _read_prepared_csv(args.annual_sediment)
        adjusted = simulate_sediment(dams, discharge, annual, start, end)
        baseline = _read_prepared_csv(args.baseline)
        daily = calculate_net_loss(baseline, adjusted)
        if _sha256(args.annual_sediment) != annual_hash:
            raise ValueError("Annual sediment input changed during calculation")
        if _sha256(args.baseline) != baseline_hash:
            raise ValueError("Baseline generation changed during calculation")
        _check_cli_inputs(args, hashes)
        metadata["model"] = "historical_hydropower_annual_sediment"
        hashes.update(annual_sediment_sha256=annual_hash, baseline_daily_sha256=baseline_hash)
        metadata["signed_net_loss"] = True
        # The paired run fingerprints both the shared operation and sediment code.
        shared_hash = bytes.fromhex(metadata["implementation_sha256"])
        sediment_hash = bytes.fromhex(_sha256(Path(__file__)))
        metadata["implementation_sha256"] = hashlib.sha256(shared_hash + sediment_hash).hexdigest()
        _write_outputs(args.output, daily, metadata)
    except (ValueError, ArithmeticError, OSError, KeyError) as exc:
        parser.exit(2, f"Input or calculation error: {exc}\n")


if __name__ == "__main__":
    sediment_cli()
