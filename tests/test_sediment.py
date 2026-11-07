"""Hand-checkable artificial tests; no study data are included.

Author: Rui Deng
"""

from __future__ import annotations

import importlib
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


def _module():
    return importlib.import_module("hma_models.sediment")


def _forcing(ids=("SYN_A",), start=1980, end=1980, flow=1.0, ssy=0.0):
    dates = pd.date_range(f"{start}-01-01", f"{end}-12-31")
    discharge = pd.DataFrame(
        [(date, dam_id, float(flow)) for date in dates for dam_id in ids],
        columns=["date", "dam_id", "flow_m3s"],
    )
    annual = pd.DataFrame(
        [(year, dam_id, float(ssy)) for year in range(start, end + 1) for dam_id in ids],
        columns=["year", "dam_id", "ssy_t_km2_yr"],
    )
    return discharge, annual


def test_first_day_uses_dead_pool_and_prior_window_and_keeps_water_sink():
    # Catches a full-pool/reset error, a current-day TE window, or removal of the sink.
    model = _module()
    dams = [model.SedimentDam("SYN_A", 1970, 100.0, 10.0, 1.0)]
    discharge, annual = _forcing(flow=1.0 / 86400, ssy=366.0)
    daily = model.simulate_historical_sediment(dams, discharge, annual, end_year=1980)
    first = daily.iloc[0]
    # Prior S=10 m3, prior outflow=1 m3/day, daily mass=1 tonne.
    expected_te = 10.0 / (0.012 * 365.25 + 1.02 * 10.0)
    assert first.trapping_efficiency == pytest.approx(expected_te)
    assert first.trapped_t == pytest.approx(expected_te)
    assert first.initial_water_added_m3 == 10.0
    assert first.water_sink_m3 == pytest.approx(expected_te / 1.2)
    assert first.stored_water_m3 == pytest.approx(11.0 - expected_te / 1.2)
    assert first.remaining_capacity_m3 == pytest.approx(100.0 - expected_te / 1.2)
    assert first.water_balance_residual_m3 == pytest.approx(0.0, abs=1e-12)
    assert first.sediment_balance_residual_t == pytest.approx(0.0, abs=1e-12)
    assert len(daily) == 366


def test_direct_upstream_sediment_routes_same_day_but_water_does_not():
    # A non-operating upstream node exports local mass; downstream gross already includes it.
    model = _module()
    dams = [
        model.SedimentDam("SYN_DOWN", 1970, 100.0, 10.0, 2.0, ("SYN_UP",)),
        model.SedimentDam("SYN_UP", 1990, 100.0, 10.0, 1.0),
    ]
    discharge, annual = _forcing(("SYN_DOWN", "SYN_UP"), flow=1 / 86400, ssy=366)
    daily = model.simulate_historical_sediment(dams, discharge, annual, end_year=1980)
    first = daily[daily.date == pd.Timestamp("1980-01-01")].set_index("dam_id")
    assert first.loc["SYN_UP", "outgoing_sediment_t"] == pytest.approx(1)
    assert first.loc["SYN_DOWN", "local_sediment_t"] == pytest.approx(1)
    assert first.loc["SYN_DOWN", "incoming_sediment_t"] == pytest.approx(2)
    assert first.loc["SYN_DOWN", "water_in_m3"] == pytest.approx(1)
    assert not bool(first.loc["SYN_UP", "is_operating"])
    assert np.isnan(first.loc["SYN_UP", "trapped_t"])
    inactive = model.summarize_annual(daily).set_index("dam_id").loc["SYN_UP"]
    assert np.isnan(inactive.trapped_t)


def test_q_squared_daily_allocation_and_negative_local_mass_clipping():
    model = _module()
    dams = [model.SedimentDam("SYN_A", 1990, 10, 0, 1)]
    discharge, annual = _forcing(flow=0, ssy=5)
    discharge.loc[0, "flow_m3s"] = 1
    discharge.loc[1, "flow_m3s"] = 2
    daily = model.simulate_historical_sediment(dams, discharge, annual, end_year=1980)
    assert daily.local_sediment_t.iloc[:3].tolist() == pytest.approx([1, 4, 0])
    assert daily.outgoing_sediment_t.sum() == pytest.approx(5)

    # Gross downstream 2 t is below the upstream 5 t, so local mass is zero.
    dams.append(model.SedimentDam("SYN_B", 1990, 10, 0, 1, ("SYN_A",)))
    discharge, annual = _forcing(("SYN_A", "SYN_B"), ssy=5)
    annual.loc[annual.dam_id == "SYN_B", "ssy_t_km2_yr"] = 2
    daily = model.simulate_historical_sediment(dams, discharge, annual, end_year=1980)
    downstream = daily[daily.dam_id == "SYN_B"]
    assert downstream.local_sediment_t.sum() == 0
    assert downstream.outgoing_sediment_t.sum() == pytest.approx(5)
    assert downstream.negative_increment_clipped.all()


def test_zero_flow_positive_annual_mass_uses_disclosed_uniform_fallback():
    model = _module()
    dams = [model.SedimentDam("SYN_A", 1990, 10, 0, 1)]
    discharge, annual = _forcing(flow=0, ssy=366)
    daily = model.simulate_historical_sediment(dams, discharge, annual, end_year=1980)
    assert daily.local_sediment_t.to_numpy() == pytest.approx(np.ones(366))
    assert daily.uniform_mass_fallback.all()


def test_capacity_bound_and_old_suspension_do_not_create_extra_deposition():
    model = _module()
    dams = [model.SedimentDam("SYN_A", 1970, 10, 10, 1)]
    discharge, annual = _forcing(flow=1 / 86400, ssy=10000)
    daily = model.simulate_historical_sediment(dams, discharge, annual, end_year=1980)
    assert daily.trapped_t.sum() == pytest.approx(12)
    assert daily.remaining_capacity_m3.min() >= 0
    assert daily.iloc[1].trapping_efficiency == 0
    assert daily.iloc[1].trapped_t == 0
    assert daily.iloc[1].outgoing_sediment_t > 0
    assert daily.sediment_balance_residual_t.abs().max() < 1e-10


def test_previous_window_rolls_over_only_after_thirty_completed_days():
    # With no release, all 30 seeded outflows have left the window on day 31.
    model = _module()
    dams = [model.SedimentDam("SYN_A", 1970, 1000, 10, 1)]
    discharge, annual = _forcing(flow=1 / 86400, ssy=0)
    daily = model.simulate_historical_sediment(dams, discharge, annual, end_year=1980)
    assert daily.residence_time_years.iloc[1] == pytest.approx(301 / (29 * 365.25))
    assert daily.residence_time_years.iloc[29] < 1e9
    assert daily.residence_time_years.iloc[30] == 1e9


def test_old_suspended_mass_does_not_settle_without_new_incoming_mass():
    model = _module()
    dams = [model.SedimentDam("SYN_A", 1970, 100, 10, 1)]
    discharge, annual = _forcing(flow=0, ssy=1)
    discharge.loc[0, "flow_m3s"] = 1 / 86400
    daily = model.simulate_historical_sediment(dams, discharge, annual, end_year=1980)
    assert daily.suspended_t.iloc[0] > 0
    assert daily.incoming_sediment_t.iloc[1] == 0
    assert daily.trapped_t.iloc[1] == 0
    assert daily.suspended_t.iloc[1] == pytest.approx(daily.suspended_t.iloc[0])


def test_sink_floor_conflict_is_exposed_in_ledger_not_reported_as_closed():
    # Positive mass / zero Q is a legacy fallback, not a physical transport claim.
    model = _module()
    dams = [model.SedimentDam("SYN_A", 1970, 100, 10, 1)]
    discharge, annual = _forcing(flow=0, ssy=3660)
    daily = model.simulate_historical_sediment(dams, discharge, annual, end_year=1980)
    assert daily.water_balance_residual_m3.iloc[0] == pytest.approx(0, abs=1e-12)
    assert daily.water_sink_m3.iloc[1] > daily.stored_water_m3.iloc[0]
    assert daily.water_balance_residual_m3.iloc[1] < -1
    assert daily.stored_water_m3.iloc[1] == 0


def test_commissioning_and_annual_cumulative_are_study_period_only():
    model = _module()
    dams = [model.SedimentDam("SYN_A", 1981, 100, 10, 1)]
    discharge, annual = _forcing(start=1980, end=1981, flow=1 / 86400, ssy=100)
    daily = model.simulate_historical_sediment(dams, discharge, annual, end_year=1981)
    assert daily[daily.year == 1980].trapped_t.isna().all()
    first = daily[daily.year == 1981].iloc[0]
    assert first.initial_water_added_m3 == 10
    annual_out = model.summarize_annual(daily).set_index("year")
    assert np.isnan(annual_out.loc[1980, "trapped_t"])
    assert annual_out.loc[1981, "cumulative_trapped_t"] == pytest.approx(
        annual_out.loc[1981, "trapped_t"]
    )

    old_dam = [model.SedimentDam("SYN_A", 1960, 100, 10, 1)]
    daily = model.simulate_historical_sediment(old_dam, discharge, annual, end_year=1981)
    annual_out = model.summarize_annual(daily)
    assert annual_out.cumulative_trapped_t.iloc[-1] == pytest.approx(annual_out.trapped_t.sum())
    assert daily.initial_water_added_m3.sum() == 10


def test_new_dam_window_is_seeded_from_first_study_year_not_its_build_year():
    model = _module()
    dams = [model.SedimentDam("SYN_A", 1981, 100, 10, 1)]
    discharge, annual = _forcing(start=1980, end=1981, flow=1 / 86400, ssy=1)
    discharge.loc[discharge.date.dt.year == 1981, "flow_m3s"] = 2 / 86400
    daily = model.simulate_historical_sediment(dams, discharge, annual, end_year=1981)
    first = daily[daily.year == 1981].iloc[0]
    assert first.residence_time_years == pytest.approx(10 / 365.25)


@pytest.mark.parametrize("bad", ["missing_day", "duplicate_day", "negative_q", "nan_q", "missing_ssy", "negative_ssy", "extra_id", "fractional_year"])
def test_incomplete_or_invalid_prepared_forcing_is_rejected(bad):
    model = _module()
    dams = [model.SedimentDam("SYN_A", 1970, 10, 1, 1)]
    discharge, annual = _forcing()
    if bad == "missing_day":
        discharge = discharge.iloc[1:]
    elif bad == "duplicate_day":
        discharge = pd.concat([discharge, discharge.iloc[:1]], ignore_index=True)
    elif bad == "negative_q":
        discharge.loc[0, "flow_m3s"] = -1
    elif bad == "nan_q":
        discharge.loc[0, "flow_m3s"] = np.nan
    elif bad == "missing_ssy":
        annual = annual.iloc[:0]
    elif bad == "negative_ssy":
        annual.loc[0, "ssy_t_km2_yr"] = -1
    elif bad == "extra_id":
        discharge.loc[0, "dam_id"] = "SYN_UNKNOWN"
    else:
        annual["year"] = 1980.5
    with pytest.raises(ValueError):
        model.simulate_historical_sediment(dams, discharge, annual, end_year=1980)


def test_invalid_topology_parameters_and_period_are_rejected():
    model = _module()
    discharge, annual = _forcing()
    good = model.SedimentDam("SYN_A", 1970, 10, 1, 1)
    for dams in (
        [good, good],
        [model.SedimentDam("SYN_A", 1970, 10, 1, 1, ("SYN_MISSING",))],
        [model.SedimentDam("SYN_A", 1970, 10, 1, 1, ("SYN_A",))],
        [model.SedimentDam("SYN_A", 1970, 0, 0, 1)],
        [model.SedimentDam("SYN_A", 1970, 10, 11, 1)],
        [model.SedimentDam("SYN_A", 1970.5, 10, 1, 1)],
        [model.SedimentDam("SYN_A", 1970, np.inf, 1, 1)],
    ):
        with pytest.raises(ValueError):
            model.simulate_historical_sediment(dams, discharge, annual, end_year=1980)
    cyclic = [
        model.SedimentDam("SYN_A", 1970, 10, 1, 1, ("SYN_B",)),
        model.SedimentDam("SYN_B", 1970, 10, 1, 1, ("SYN_A",)),
    ]
    with pytest.raises(ValueError, match="cycle"):
        model.simulate_historical_sediment(cyclic, discharge, annual, end_year=1980)
    for start, end in ((1979, 1980), (1980, 2020), (1981, 1980)):
        with pytest.raises(ValueError):
            model.simulate_historical_sediment([good], discharge, annual, start_year=start, end_year=end)


@pytest.mark.parametrize("opaque_id", ["SYN_A", "NA", "0007"])
def test_cli_writes_artificial_outputs_and_refuses_nonempty_destination(tmp_path, opaque_id):
    # Catches broken prepared-CSV parsing and accidental output overwrites.
    dams = pd.DataFrame([{
        "dam_id": opaque_id, "built_year": 1970, "total_capacity_m3": 100,
        "dead_capacity_m3": 10, "basin_area_km2": 1, "upstream_ids": "",
    }])
    discharge, annual = _forcing(ids=(opaque_id,), flow=1 / 86400, ssy=1)
    for name, table in (("dams", dams), ("discharge", discharge), ("ssy", annual)):
        table.to_csv(tmp_path / f"{name}.csv", index=False)
    script = Path(__file__).resolve().parents[1] / "scripts" / "02_historical_sediment.py"
    command = [sys.executable, str(script), "--dams", str(tmp_path / "dams.csv"),
               "--discharge", str(tmp_path / "discharge.csv"), "--ssy", str(tmp_path / "ssy.csv"),
               "--output", str(tmp_path / "out"), "--start-year", "1980", "--end-year", "1980"]
    run = subprocess.run(command, capture_output=True, text=True)
    assert run.returncode == 0, run.stderr
    assert len(pd.read_csv(tmp_path / "out" / "daily_sediment.csv")) == 366
    result = pd.read_csv(tmp_path / "out" / "annual_sediment.csv", dtype={"dam_id": str}, keep_default_na=False)
    assert result.dam_id.tolist() == [opaque_id]
    original = (tmp_path / "out" / "annual_sediment.csv").read_bytes()
    rerun = subprocess.run(command, capture_output=True, text=True)
    assert rerun.returncode != 0
    assert (tmp_path / "out" / "annual_sediment.csv").read_bytes() == original


@pytest.mark.parametrize("target", ["dams", "discharge", "ssy"])
@pytest.mark.parametrize("bad_header", ["duplicate", "empty", "whitespace"])
def test_cli_rejects_raw_duplicate_or_empty_csv_headers_before_pandas_mangling(
    tmp_path, target, bad_header,
):
    # Catches silently selecting the first duplicate column after pandas renames it.
    import csv

    dams = pd.DataFrame([{
        "dam_id": "SYN_A", "built_year": 1970, "total_capacity_m3": 100,
        "dead_capacity_m3": 10, "basin_area_km2": 1, "upstream_ids": "",
    }])
    discharge, annual = _forcing(flow=1 / 86400, ssy=1)
    tables = {"dams": dams, "discharge": discharge, "ssy": annual}
    duplicate_names = {"dams": "total_capacity_m3", "discharge": "flow_m3s", "ssy": "ssy_t_km2_yr"}
    for name, table in tables.items():
        path = tmp_path / f"{name}.csv"
        if name != target:
            table.to_csv(path, index=False)
            continue
        extra = duplicate_names[name] if bad_header == "duplicate" else ("" if bad_header == "empty" else "   ")
        with path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow([*table.columns, extra])
            writer.writerows([*row, 999] for row in table.itertuples(index=False, name=None))
    script = Path(__file__).resolve().parents[1] / "scripts" / "02_historical_sediment.py"
    command = [sys.executable, str(script), "--dams", str(tmp_path / "dams.csv"),
               "--discharge", str(tmp_path / "discharge.csv"), "--ssy", str(tmp_path / "ssy.csv"),
               "--output", str(tmp_path / "out"), "--end-year", "1980"]
    result = subprocess.run(command, capture_output=True, text=True)
    assert result.returncode != 0, "Malformed CSV headers must not be silently accepted"
    assert "header" in result.stderr.lower() or "column" in result.stderr.lower()
    assert not (tmp_path / "out").exists()
