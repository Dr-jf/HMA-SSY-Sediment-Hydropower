# Sediment yield, reservoir trapping and hydropower in High Mountain Asia

This repository contains four Python examples for specific sediment yield (SSY), reservoir sediment trapping and historical hydropower generation in High Mountain Asia, together with the supplementary inventory of 325 reservoirs and 179 hydropower stations.

## Scripts

| Script | Calculation |
| --- | --- |
| [`ssy_model.py`](ssy_model.py) | Fit the 17-predictor XGBoost ensemble on square-root SSY. The default 70/30 split of 1,963 records gives 1,374 training and 589 validation records. |
| [`reservoir_sediment.py`](reservoir_sediment.py) | Calculate daily and annual trapped sediment for the 325-reservoir network. |
| [`hydropower_baseline.py`](hydropower_baseline.py) | Calculate daily generation with fixed original active storage. |
| [`hydropower_sedimentation.py`](hydropower_sedimentation.py) | Calculate daily generation with annual active-storage updates and pair it with the baseline to estimate generation loss. |

## Reservoir inventory

[Supplementary Data - HMA reservoir inventory.xlsx](Supplementary%20Data%20-%20HMA%20reservoir%20inventory.xlsx) contains 325 reservoir records, including 179 in the `Hydropower & irrigation` category and 146 in the `Irrigation` category. The `Reservoir inventory` worksheet holds the data, and the workbook's `README` worksheet defines the 19 fields.

| Field | Parameter or description | Unit | Records with values |
| --- | --- | --- | --- |
| `damkey` | Unique study identifier | — | 325 |
| `dam_name` | English reservoir or dam label used in the model inventory | — | 325 |
| `longitude_deg_e` | Longitude of the registered model point | °E | 325 |
| `latitude_deg_n` | Latitude of the registered model point | °N | 325 |
| `construction_year` | Recorded dam construction year | year | 325 |
| `mainuse` | Analysis category | — | 325 |
| `reservoir_area_km2` | Reservoir surface area used in the model | km² | 325 |
| `total_storage_km3` | Initial total storage, equal to active plus dead storage | km³ | 325 |
| `active_storage_km3` | Initial active storage above the dead-water level | km³ | 325 |
| `dead_storage_km3` | Initial dead storage below the dead-water level | km³ | 325 |
| `normal_operating_level_m` | Normal pool elevation, the upper operating water level | m | 325 |
| `dead_water_level_m` | Dead-water level, the lower operating water level | m | 325 |
| `morphology_m` | Morphology exponent for the volume–depth relation | dimensionless | 325 |
| `max_generation_head_m` | Maximum modelled generating head at the normal operating level | m | 179 |
| `installed_capacity_mw` | Installed electrical generating capacity used in the model | MW | 179 |
| `max_turbine_flow_m3_s` | Maximum turbine discharge imposed in the hydropower simulation | m³ s⁻¹ | 179 |
| `annual_generation_ref_gwh_yr` | Reported annual generation reference | GWh yr⁻¹ | 158 |
| `annual_sediment_retention_ref_mt_yr` | Reported or derived comparable annual sediment mass retained by the reservoir | Mt yr⁻¹ | 53 |
| `reference_url` | Source reference URL | — | 325 |

Water-level elevations retain the source site's datum. The morphology exponent describes the volume–depth relation and uses the geometric estimate of surface area × maximum depth proxy / initial total volume, with all quantities in SI units. A slash (`/`) marks an unavailable or inapplicable value.

To prepare the model CSVs described below, map `damkey` to `dam_id` and convert storage from km³ to m³ by multiplying by 10⁹. Use the recorded construction year consistently with the scripts' activation-year convention. Prepare the remaining model columns, reservoir-network relationships and forcing series using the field definitions below.

## Running the examples

Use Python 3.11 or 3.12. Install the dependencies and run the scripts from the repository directory.

```sh
python -m pip install numpy==1.26.4 pandas==2.2.3 scipy==1.12.0 scikit-learn==1.5.2 xgboost==2.1.3
python ssy_model.py --input inputs/ssy.csv --output outputs/ssy
python reservoir_sediment.py --dams inputs/sediment_dams.csv --discharge inputs/discharge.csv --ssy inputs/annual_ssy.csv --output outputs/sediment
python hydropower_baseline.py --dams inputs/hydropower_dams.csv --discharge inputs/discharge.csv --output outputs/baseline
python hydropower_sedimentation.py --dams inputs/hydropower_dams.csv --discharge inputs/discharge.csv --annual-sediment outputs/sediment/annual_sediment.csv --baseline outputs/baseline/daily_generation.csv --output outputs/sedimentation
```

Use a new or empty directory for each output. Keep the baseline's `manifest.json` next to its `daily_generation.csv` so the paired run can check the inputs and dates. Keep both hydropower scripts in the same directory.

## Input tables

Supply UTF-8 CSVs with one header row. Keep `dam_id` consistent across tables, including leading zeros.

- **SSY training.** One row per record with `SSY` (t km⁻² yr⁻¹) and the predictors `PGA`, `S`, `V`, `P`, `T`, `L`, `Ksn`, `A`, `G`, `F`, `R`, `GMC`, `Q`, `SC`, `SM`, `DeltaG` and `DeltaSC`. Use the study's predictor definitions and lithology encoding. The annual runoff predictor `Q` is in mm yr⁻¹. Missing predictors may be left blank. `SSY` must be finite and non-negative.
- **Sediment inventory.** `dam_id`, `built_year`, `total_capacity_m3`, `dead_capacity_m3`, `basin_area_km2` and `upstream_ids`. List direct upstream reservoir IDs separated by semicolons. Leave the field blank for headwater reservoirs.
- **Daily discharge.** `date` (YYYY-MM-DD), `dam_id` and `flow_m3s` (m³ s⁻¹). Supply one value per required dam and day.
- **Annual SSY for trapping.** `year`, `dam_id` and `ssy_t_km2_yr` for the full upstream drainage area of each reservoir. Prepare these basin means from the spatial SSY inputs before running the trapping example.
- **Hydropower inventory.** `dam_id`, `built_year`, `active_capacity_m3`, `dead_capacity_m3`, `reservoir_area_km2`, `dam_height_m`, `normal_level_m`, `dead_level_m`, `intake_level_m`, `installed_capacity_mw`, `max_discharge_m3s` and `wet_months`. List consecutive wet-season months within the calendar year separated by semicolons. The optional `efficiency` defaults to 0.90.

The historical examples start from original storage capacity in 1980 or the commissioning year, whichever is later. The trapping example allocates incremental annual sediment using daily discharge squared and routes sediment through the supplied reservoir network. It uses a 30-day window for the Brune–Gill M trapping curve and deducts deposited volume from stored water before release.

The hydropower examples use identical inventories and discharge. The sediment run takes annual increments from `trapped_t` in `annual_sediment.csv`. It applies each calendar year's accumulated sediment after generation on the last day of the wet season. The new capacity takes effect the following day. Both scripts report daily and annual generation in MWh. The sediment run reports signed loss as baseline generation minus sediment-scenario generation.

## License

MIT License

Copyright (c) 2026 Rui Deng

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
