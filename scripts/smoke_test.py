"""Exercise four entry points using temporary, wholly artificial inputs.

No research data are read, shipped or inferred. Author: Rui Deng.
"""

from pathlib import Path
import json
import os
import subprocess
import sys
import tempfile
import time

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from hma_models.ssy import FEATURES, predict_ssy


def run_entry(name, *arguments):
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / name), *map(str, arguments)],
        cwd=ROOT, env=environment, capture_output=True, text=True, timeout=180,
    )
    if result.returncode:
        raise RuntimeError(f"Artificial smoke run failed: {name}\n{result.stdout}\n{result.stderr}")
    print(f"PASS {name}", flush=True)


def main():
    started = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix="hma_artificial_smoke_") as temporary:
        folder = Path(temporary)
        rng = np.random.default_rng(1701)
        station = pd.DataFrame(rng.normal(size=(40, len(FEATURES))), columns=FEATURES)
        station["SSY"] = (10.0 + 2.0 * np.maximum(station["P"], 0.0)) ** 2
        station_path = folder / "artificial_ssy.csv"
        station.to_csv(station_path, index=False)
        model_dir = folder / "artificial_model"
        run_entry("01_train_ssy.py", "--input", station_path, "--output", model_dir,
                  "--expected-samples", "40", "--threads", "1")
        metadata = json.loads((model_dir / "metadata.json").read_text(encoding="utf-8"))
        assert metadata["sample_counts"] == {"total": 40, "training": 28, "test": 12}
        prediction = predict_ssy(model_dir, station)
        assert prediction.shape == (40,) and np.isfinite(prediction).all() and (prediction >= 0).all()

        # Literal fixture parameters are artificial, not reservoir measurements.
        dates = pd.date_range("1980-01-01", "1980-12-31", freq="D")
        flow_path = folder / "artificial_discharge.csv"
        pd.DataFrame({
            "date": dates.strftime("%Y-%m-%d"), "dam_id": "1",
            "flow_m3s": 2.0 + np.sin(np.arange(len(dates)) * 2 * np.pi / len(dates)),
        }).to_csv(flow_path, index=False)
        sediment_dams = folder / "artificial_sediment_dams.csv"
        pd.DataFrame([dict(
            dam_id="1", built_year=1970, total_capacity_m3=11000.0,
            dead_capacity_m3=1000.0, basin_area_km2=1.0, upstream_ids="",
        )]).to_csv(sediment_dams, index=False)
        annual_ssy = folder / "artificial_annual_ssy.csv"
        pd.DataFrame([dict(year=1980, dam_id="1", ssy_t_km2_yr=100000.0)]).to_csv(annual_ssy, index=False)
        sediment_output = folder / "artificial_sediment_output"
        run_entry("02_historical_sediment.py", "--dams", sediment_dams,
                  "--discharge", flow_path, "--ssy", annual_ssy, "--output", sediment_output,
                  "--start-year", "1980", "--end-year", "1980")
        annual_mass = pd.read_csv(sediment_output / "annual_sediment.csv")
        assert len(annual_mass) == 1 and annual_mass.loc[0, "trapped_t"] > 0

        hydropower_dams = folder / "artificial_hydropower_dams.csv"
        pd.DataFrame([dict(
            dam_id="1", built_year=1970, active_capacity_m3=10000.0,
            dead_capacity_m3=1000.0, reservoir_area_km2=0.05, dam_height_m=20.0,
            normal_level_m=120.0, dead_level_m=110.0, intake_level_m=100.0,
            installed_capacity_mw=0.5, max_discharge_m3s=2.5,
            wet_months="6;7;8;9", efficiency=0.9,
        )]).to_csv(hydropower_dams, index=False)
        baseline_output = folder / "artificial_baseline"
        adjusted_output = folder / "artificial_adjusted"
        common = ("--dams", hydropower_dams, "--discharge", flow_path,
                  "--start-year", "1980", "--end-year", "1980")
        run_entry("03_historical_hydropower_baseline.py", *common, "--output", baseline_output)
        run_entry("04_historical_hydropower_sediment.py", *common,
                  "--annual-sediment", sediment_output / "annual_sediment.csv",
                  "--baseline", baseline_output / "daily_generation.csv", "--output", adjusted_output)
        baseline = pd.read_csv(baseline_output / "daily_generation.csv")
        adjusted = pd.read_csv(adjusted_output / "daily_generation.csv")
        assert len(baseline) == len(adjusted) == 366
        assert np.array_equal(baseline["date"], adjusted["date"])
        assert np.allclose(adjusted["net_loss_mwh"], baseline["energy_mwh"] - adjusted["energy_mwh"])
        assert adjusted["active_capacity_end_m3"].min() < 10000.0
        assert adjusted["water_balance_residual_m3"].abs().max() <= 1e-5

    print(f"Artificial smoke checks completed in {time.perf_counter() - started:.2f} s; no manuscript results calculated.")


if __name__ == "__main__":
    main()
