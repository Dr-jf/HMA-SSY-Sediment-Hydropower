"""Artificial water, energy, annual-update and counterfactual checks."""

import importlib
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest


def api():
    return importlib.import_module("hma_models.hydropower")


def dam(**changes):
    values = dict(
        dam_id="artificial-a", built_year=1970, active_capacity_m3=1000.0,
        dead_capacity_m3=1000.0, reservoir_area_km2=0.0001,
        dam_height_m=20.0, normal_level_m=120.0, dead_level_m=110.0,
        intake_level_m=100.0, installed_capacity_mw=100.0,
        max_discharge_m3s=1.0, wet_months=(6, 7, 8, 9), efficiency=0.9,
    )
    values.update(changes)
    return api().HydropowerDam(**values)


def forcing(start="1980-01-01", end="1980-12-31", flow=0.2, dam_id="artificial-a"):
    dates = pd.date_range(start, end, freq="D")
    return pd.DataFrame({"date": dates, "dam_id": dam_id, "flow_m3s": float(flow)})


def annual(mass=0.0):
    return pd.DataFrame({"dam_id": ["artificial-a"], "year": [1980], "trapped_t": [mass]})


def assert_water_balance(result):
    independent = (
        result.storage_end_m3 - result.storage_start_m3
        - 86400.0 * result.flow_m3s + 86400.0 * result.routed_m3s
        + result.spill_m3
    )
    assert np.abs(independent).max() < 1e-7
    assert (result.storage_end_m3 >= 0).all()
    assert (result.storage_end_m3 <= result.active_capacity_end_m3).all()


def test_zero_sediment_gives_identical_daily_counterfactual_and_calendar():
    baseline = api().simulate_baseline([dam()], forcing(), 1980, 1980)
    sediment = api().simulate_sediment([dam()], forcing(), annual(), 1980, 1980)
    assert len(baseline) == 366
    assert (baseline.date == pd.Timestamp("1980-02-29")).sum() == 1
    np.testing.assert_array_equal(baseline.energy_mwh, sediment.energy_mwh)
    assert baseline.iloc[0].storage_start_m3 == 1000
    assert_water_balance(baseline)
    assert_water_balance(sediment)


def test_dry_release_uses_remaining_days_and_analytical_linear_head():
    reservoir = dam(active_capacity_m3=2000, dead_capacity_m3=2000, reservoir_area_km2=0.0002)
    result = api().simulate_baseline([reservoir], forcing(flow=0), 1980, 1980)
    first = result.iloc[0]
    # Leap-year January through May has 152 days; m=1 for this geometry.
    expected_q = 2000 / (152 * 86400)
    expected_energy = 24 * 0.00981 * 0.9 * expected_q * (20 - 5 / 152)
    assert first.turbine_m3s == pytest.approx(expected_q)
    assert first.storage_end_m3 == pytest.approx(2000 * 151 / 152)
    assert first.energy_mwh == pytest.approx(expected_energy)
    assert_water_balance(result)


def test_power_is_capped_during_head_integration():
    reservoir = dam(installed_capacity_mw=1.3, max_discharge_m3s=10)
    actual, head = api().daily_energy_mwh(reservoir, 1000, 1000, 0, 10)
    # Independently integrate the clipped linear-head power curve.
    heads = np.linspace(10, 20, 100001)
    power = np.minimum(1.3, 0.00981 * 0.9 * 10 * heads)
    expected = 24 * np.trapz(power, heads) / 10
    assert actual == pytest.approx(expected, abs=1e-7)
    assert head == pytest.approx(15)


def test_installed_power_and_maximum_discharge_hold_with_excess_inflow():
    result = api().simulate_baseline(
        [dam(installed_capacity_mw=0.01)], forcing(flow=10), 1980, 1980
    )
    assert result.energy_mwh.max() == pytest.approx(0.24)
    assert result.turbine_m3s.max() <= 1
    assert result.spill_m3.min() == pytest.approx(9 * 86400)
    assert_water_balance(result)


def test_wet_refill_keeps_limited_inflow_and_allows_unfilled_storage():
    q = forcing(flow=0)
    q.loc[q.date.dt.month.isin([6, 7, 8, 9]), "flow_m3s"] = 1e-6
    result = api().simulate_baseline([dam()], q, 1980, 1980)
    wet = result[result.date.dt.month.isin([6, 7, 8, 9])]
    assert wet.energy_mwh.sum() == 0
    assert wet.iloc[-1].storage_end_m3 == pytest.approx(122 * 86400 * 1e-6, abs=1e-7)
    assert_water_balance(result)


def test_whole_year_mass_is_applied_after_wet_end_generation_and_spills_water():
    q = forcing(flow=10)
    baseline = api().simulate_baseline([dam()], q, 1980, 1980)
    result = api().simulate_sediment([dam()], q, annual(1200), 1980, 1980)
    before = result[result.date == "1980-09-29"].iloc[0]
    update = result[result.date == "1980-09-30"].iloc[0]
    after = result[result.date == "1980-10-01"].iloc[0]
    assert before.active_capacity_end_m3 == 1000
    assert update.active_capacity_start_m3 == 1000
    assert update.active_capacity_end_m3 == 500
    assert update.capacity_update_spill_m3 == 500
    assert update.energy_mwh == baseline[baseline.date == "1980-09-30"].iloc[0].energy_mwh
    assert after.active_capacity_start_m3 == 500
    assert after.storage_start_m3 == 500
    assert_water_balance(result)


def test_mass_is_accumulated_as_annual_increments_without_double_counting():
    q = forcing(end="1981-12-31", flow=10)
    masses = pd.DataFrame({"dam_id": ["artificial-a"] * 2, "year": [1980, 1981], "trapped_t": [600, 600]})
    result = api().simulate_sediment([dam()], q, masses, 1980, 1981)
    assert result.iloc[-1].active_capacity_end_m3 == 500
    assert result.iloc[-1].cumulative_sediment_m3 == 1000
    assert result.iloc[-1].cumulative_dead_loss_m3 == 500


def test_full_capacity_exhaustion_routes_inflow_with_normal_pool_head():
    q = forcing(flow=0.2)
    result = api().simulate_sediment([dam()], q, annual(2400), 1980, 1980)
    after = result[result.date == "1980-10-01"].iloc[0]
    assert after.active_capacity_start_m3 == 0
    assert after.storage_start_m3 == after.storage_end_m3 == 0
    assert after.head_mean_m == 20
    assert after.energy_mwh == pytest.approx(24 * 0.00981 * 0.9 * 0.2 * 20)
    assert_water_balance(result)


def test_zero_dead_height_and_redistribution_do_not_create_negative_capacity():
    reservoir = dam(dead_level_m=100, intake_level_m=90)
    allocation = api().allocate_cumulative_sediment(reservoir, 1500)
    assert allocation["active_loss_m3"] == 1000
    assert allocation["dead_loss_m3"] == 500
    assert allocation["unallocated_m3"] == 0
    full = api().allocate_cumulative_sediment(reservoir, 2300)
    assert full["unallocated_m3"] == 300


def test_period_reset_has_no_prehistory_and_preconstruction_days_are_excluded():
    q = forcing(end="1981-12-31")
    both = pd.concat([q, q.assign(dam_id="artificial-b")], ignore_index=True)
    dams = [dam(), dam(dam_id="artificial-b", built_year=1981)]
    result = api().simulate_baseline(dams, both, 1980, 1981)
    new = result[result.dam_id == "artificial-b"]
    assert len(new) == 365
    assert new.iloc[0].date == pd.Timestamp("1981-01-01")
    assert new.iloc[0].storage_start_m3 == 1000
    a = api().simulate_baseline([dam(built_year=1970)], q, 1981, 1981)
    b = api().simulate_baseline([dam(built_year=1981)], q, 1981, 1981)
    np.testing.assert_array_equal(a.energy_mwh, b.energy_mwh)


def test_inactive_and_other_dam_sediment_rows_may_be_na():
    masses = pd.DataFrame({"dam_id": ["artificial-a", "other"], "year": [1980, 1980], "trapped_t": [np.nan, np.nan]})
    result = api().simulate_sediment([dam(built_year=1981)], forcing(), masses, 1980, 1980)
    assert result.empty


@pytest.mark.parametrize("months", [(1, 2, 12), (6, 8), (), tuple(range(1, 13))])
def test_invalid_annual_wet_period_is_rejected(months):
    with pytest.raises(ValueError, match="wet"):
        dam(wet_months=months)


@pytest.mark.parametrize("fault", ["missing", "duplicate", "negative", "infinite", "unknown_id"])
def test_invalid_operational_discharge_is_rejected(fault):
    q = forcing()
    if fault == "missing":
        q = q.drop(index=59)
    elif fault == "duplicate":
        q = pd.concat([q, q.iloc[[0]]], ignore_index=True)
    elif fault == "negative":
        q.loc[0, "flow_m3s"] = -1
    elif fault == "infinite":
        q.loc[0, "flow_m3s"] = np.inf
    else:
        q.loc[0, "dam_id"] = "unknown"
    with pytest.raises(ValueError):
        api().simulate_baseline([dam()], q, 1980, 1980)


@pytest.mark.parametrize("fault", ["negative", "missing", "duplicate", "nan"])
def test_invalid_operational_annual_increment_is_rejected(fault):
    mass = annual(-1 if fault == "negative" else 0)
    if fault == "missing":
        mass = mass.iloc[:0]
    elif fault == "duplicate":
        mass = pd.concat([mass, mass], ignore_index=True)
    elif fault == "nan":
        mass.loc[0, "trapped_t"] = np.nan
    with pytest.raises(ValueError):
        api().simulate_sediment([dam()], forcing(), mass, 1980, 1980)


def test_historical_window_and_unique_dam_ids_are_enforced():
    with pytest.raises(ValueError):
        api().simulate_baseline([dam()], forcing(), 1979, 1980)
    with pytest.raises(ValueError):
        api().simulate_baseline([dam()], forcing(), 1980, 2020)
    with pytest.raises(ValueError):
        api().simulate_baseline([dam(), dam()], forcing(), 1980, 1980)


def test_signed_net_loss_preserves_gains_and_rejects_unmatched_or_duplicate_keys():
    base = pd.DataFrame({"date": ["1980-01-01", "1980-01-02"], "dam_id": ["artificial-a"] * 2, "energy_mwh": [2.0, 2.0]})
    adjusted = base.assign(energy_mwh=[3.0, 1.0])
    result = api().calculate_net_loss(base, adjusted)
    assert result.net_loss_mwh.tolist() == [-1, 1]
    with pytest.raises(ValueError):
        api().calculate_net_loss(base.iloc[:1], adjusted)
    with pytest.raises(ValueError):
        api().calculate_net_loss(pd.concat([base, base.iloc[:1]]), adjusted)


def write_cli_inputs(tmp_path):
    reservoir = dam()
    # Synthetic fixture files exist only under pytest's temporary directory.
    params = {key: value for key, value in reservoir.__dict__.items() if not key.startswith("_")}
    params["wet_months"] = "6;7;8;9"
    dam_path, flow_path, sediment_path = [tmp_path / name for name in ["dams.csv", "flow.csv", "annual.csv"]]
    pd.DataFrame([params]).to_csv(dam_path, index=False)
    forcing().to_csv(flow_path, index=False, date_format="%Y-%m-%d")
    annual(1200).to_csv(sediment_path, index=False)
    return dam_path, flow_path, sediment_path


def run_cli(script_name, *arguments):
    root = Path(__file__).resolve().parents[1]
    return subprocess.run([sys.executable, str(root / "scripts" / script_name), *map(str, arguments)], cwd=root, capture_output=True, text=True)


def test_cli_pair_provenance_alignment_privacy_and_no_overwrite(tmp_path):
    dams, q, mass = write_cli_inputs(tmp_path)
    common = ["--dams", dams, "--discharge", q, "--start-year", 1980, "--end-year", 1980]
    baseline_dir, sediment_dir = tmp_path / "baseline", tmp_path / "sediment"
    base_run = run_cli("03_historical_hydropower_baseline.py", *common, "--output", baseline_dir)
    assert base_run.returncode == 0, base_run.stderr
    assert (baseline_dir / "daily_generation.csv").is_file()
    assert str(tmp_path) not in (baseline_dir / "manifest.json").read_text()
    sed_run = run_cli("04_historical_hydropower_sediment.py", *common, "--annual-sediment", mass, "--baseline", baseline_dir / "daily_generation.csv", "--output", sediment_dir)
    assert sed_run.returncode == 0, sed_run.stderr
    result = pd.read_csv(sediment_dir / "daily_generation.csv")
    assert len(result) == 366
    assert "net_loss_mwh" in result
    repeated = run_cli("03_historical_hydropower_baseline.py", *common, "--output", baseline_dir)
    assert repeated.returncode != 0
    # A different counterfactual forcing must be rejected even with identical keys.
    changed = pd.read_csv(q)
    changed.loc[0, "flow_m3s"] = 0.3
    changed.to_csv(q, index=False)
    bad = run_cli("04_historical_hydropower_sediment.py", *common, "--annual-sediment", mass, "--baseline", baseline_dir / "daily_generation.csv", "--output", tmp_path / "mismatched")
    assert bad.returncode != 0
    assert not (tmp_path / "mismatched" / "daily_generation.csv").exists()


def test_cli_rejects_tampered_or_incomplete_baseline(tmp_path):
    dams, q, mass = write_cli_inputs(tmp_path)
    common = ["--dams", dams, "--discharge", q, "--start-year", 1980, "--end-year", 1980]
    folder = tmp_path / "baseline"
    assert run_cli("03_historical_hydropower_baseline.py", *common, "--output", folder).returncode == 0
    baseline = pd.read_csv(folder / "daily_generation.csv").iloc[:-1]
    baseline.to_csv(folder / "daily_generation.csv", index=False)
    bad = run_cli("04_historical_hydropower_sediment.py", *common, "--annual-sediment", mass, "--baseline", folder / "daily_generation.csv", "--output", tmp_path / "bad")
    assert bad.returncode != 0


@pytest.mark.parametrize("role,duplicate", [(r, d) for r in ["dams", "discharge", "annual", "baseline"] for d in [True, False]])
def test_cli_rejects_raw_duplicate_and_empty_headers(tmp_path, role, duplicate):
    dams, q, mass = write_cli_inputs(tmp_path)
    common = ["--dams", dams, "--discharge", q, "--start-year", 1980, "--end-year", 1980]
    baseline_dir = tmp_path / "baseline"
    if role in {"annual", "baseline"}:
        completed = run_cli("03_historical_hydropower_baseline.py", *common, "--output", baseline_dir)
        assert completed.returncode == 0, completed.stderr
    target = {"dams": dams, "discharge": q, "annual": mass, "baseline": baseline_dir / "daily_generation.csv"}[role]
    name = {"dams": "active_capacity_m3", "discharge": "flow_m3s", "annual": "trapped_t", "baseline": "energy_mwh"}[role]
    with target.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.reader(stream))
    rows[0].append(name if duplicate else "")
    for row in rows[1:]:
        row.append("999")
    with target.open("w", newline="", encoding="utf-8") as stream:
        csv.writer(stream).writerows(rows)
    if role == "baseline":
        # Bypass the separate integrity check to exercise the parser boundary.
        manifest_path = baseline_dir / "manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["daily_csv_sha256"] = hashlib.sha256(target.read_bytes()).hexdigest()
        manifest_path.write_text(json.dumps(manifest))
    if role in {"dams", "discharge"}:
        result = run_cli("03_historical_hydropower_baseline.py", *common, "--output", tmp_path / "new")
    else:
        result = run_cli("04_historical_hydropower_sediment.py", *common, "--annual-sediment", mass, "--baseline", baseline_dir / "daily_generation.csv", "--output", tmp_path / "new")
    assert result.returncode != 0
    assert "header" in result.stderr.lower()


def test_cli_preserves_opaque_na_identifier_across_both_scenarios(tmp_path):
    dams, q, mass = write_cli_inputs(tmp_path)
    for path in [dams, q, mass]:
        frame = pd.read_csv(path, keep_default_na=False)
        frame["dam_id"] = "NA"
        frame.to_csv(path, index=False)
    increments = pd.read_csv(mass, keep_default_na=False)
    increments = pd.concat([increments, pd.DataFrame({"dam_id": ["other", "NA"], "year": [1980, 1981], "trapped_t": ["", ""]})], ignore_index=True)
    increments.to_csv(mass, index=False)
    common = ["--dams", dams, "--discharge", q, "--start-year", 1980, "--end-year", 1980]
    baseline_dir = tmp_path / "baseline"
    baseline = run_cli("03_historical_hydropower_baseline.py", *common, "--output", baseline_dir)
    assert baseline.returncode == 0, baseline.stderr
    adjusted = run_cli("04_historical_hydropower_sediment.py", *common, "--annual-sediment", mass, "--baseline", baseline_dir / "daily_generation.csv", "--output", tmp_path / "adjusted")
    assert adjusted.returncode == 0, adjusted.stderr
    saved = pd.read_csv(tmp_path / "adjusted" / "daily_generation.csv", keep_default_na=False)
    assert saved.dam_id.unique().tolist() == ["NA"]
    assert len(saved) == 366
