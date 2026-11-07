# Private input schemas and units

All tables below are **supplied locally by the person running the code**. None is distributed in this repository. CSV is the interchange format; use UTF-8 and one header row. Generated outputs should remain private as well.

## 1. SSY training

The training CSV contains one row per observation, 17 numeric predictor columns in this order, and a numeric `SSY` target:

| Column | Meaning / required representation |
| --- | --- |
| `PGA` | Peak ground acceleration predictor, in the original source-product units |
| `S` | Mean slope, degrees |
| `V` | NDVI, dimensionless; not unconverted LAI |
| `P` | Annual precipitation, mm per year |
| `T` | Annual air temperature, degrees Celsius |
| `L` | Lithology predictor, in the original numerical encoding |
| `Ksn` | Normalized channel-steepness predictor, using the original source definition and units |
| `A` | Drainage-basin area, km2 |
| `G` | Glacier cover, percent |
| `F` | Permafrost cover, percent |
| `R` | Relief, m |
| `GMC` | Glacier mass-change predictor, m water equivalent per year |
| `Q` | Annual runoff depth, mm per year; not daily streamflow in m3/s |
| `SC` | Snow cover, percent |
| `SM` | Snowmelt, mm water equivalent per year |
| `DeltaG` | Annual glacier-cover change, percentage points per year |
| `DeltaSC` | Annual snow-cover change, percentage points per year |
| `SSY` | Specific sediment yield, tonnes per km2 per year |

Predictor missing values may be represented as empty CSV cells and are handled by XGBoost's native missing-value routing. Infinite or non-numeric predictors are invalid. Targets must be finite and non-negative; missing targets are not imputed. `DeltaG`/`DeltaSC` are the ASCII interface names; supported Greek aliases are documented in `ssy.md`.

The default training command validates exactly 1,963 rows. The code does not reconstruct the extra 52 records, infer station identities, or verify the scientific provenance of a private table. Preserve the source definitions and preprocessing used to construct the predictors. This core release does not extract environmental rasters or resolve ambiguities in source-product units.

## 2. Historical sediment inputs

### Reservoir attributes

| Column | Unit / requirement |
| --- | --- |
| `dam_id` | Unique opaque string key within the supplied tables; retain leading zeros consistently |
| `built_year` | Commissioning year; activation is 1 January |
| `total_capacity_m3` | Original total capacity, m3, positive |
| `dead_capacity_m3` | Original dead capacity, m3, between zero and total capacity |
| `basin_area_km2` | Complete upstream controlling-basin area, km2, positive |
| `upstream_ids` | Semicolon-separated keys of direct upstream reservoirs; empty for no upstream reservoir |

Provide direct relationships, not every upstream ancestor. IDs must exist in the supplied inventory, be unique, and form an acyclic network. Reservoirs that are not yet commissioned may still pass sediment through the network; their unoperated trapping metrics are not observations of zero.

### Daily discharge

| Column | Unit / requirement |
| --- | --- |
| `date` | Calendar date, `YYYY-MM-DD` |
| `dam_id` | Matching reservoir key |
| `flow_m3s` | Direct nearest-reach discharge, m3/s, finite and non-negative |

Provide a complete, duplicate-free daily series for the required historical dates. Each reservoir uses its own direct discharge. Do not add upstream model outflow to an already complete downstream forcing series.

### Annual SSY

| Column | Unit / requirement |
| --- | --- |
| `year` | Calendar year |
| `dam_id` | Matching reservoir key |
| `ssy_t_km2_yr` | Complete-control-basin mean SSY, tonnes per km2 per year, finite and non-negative |

These are already prepared basin means. The code multiplies by complete controlling-basin area and then subtracts direct-upstream gross annual sediment to obtain the local increment. Supplying an already incremental SSY product would subtract upstream sediment twice.

## 3. Hydropower attributes

| Column | Unit / requirement |
| --- | --- |
| `dam_id` | Unique opaque string key, encoded consistently across inputs |
| `built_year` | Commissioning year |
| `active_capacity_m3` | Original active capacity, m3, positive |
| `dead_capacity_m3` | Original dead capacity, m3, non-negative |
| `reservoir_area_km2` | Reservoir area, km2, positive |
| `dam_height_m` | Dam height, m, positive |
| `normal_level_m` | Normal-pool elevation, m |
| `dead_level_m` | Dead-pool elevation, m, below normal level |
| `intake_level_m` | Intake elevation, m; operating value is limited to the dead-pool level |
| `installed_capacity_mw` | Installed capacity, MW, positive |
| `max_discharge_m3s` | Maximum turbine discharge, m3/s, positive |
| `wet_months` | Semicolon-separated month numbers defining one continuous, non-cross-year wet season |
| `efficiency` | Dimensionless; default 0.90 if omitted |

The discharge table has the same `date,dam_id,flow_m3s` schema as above. The active hydropower inventory may be a subset of the sediment inventory. Names and locations are not required public interfaces.

## 4. Sediment input to the hydropower update

Use the historical sediment command's `annual_sediment.csv`, with these fields:

| Column | Unit / requirement |
| --- | --- |
| `dam_id` | Matching hydropower key |
| `year` | Calendar year |
| `trapped_t` | **Annual newly trapped mass**, tonnes; not cumulative mass, Mt, volume or a percentage |

Other sediment-output columns are not required by the hydropower model. The sediment table can contain the full inventory, including uncommissioned entries, but each selected hydropower reservoir needs valid annual increments for its operating years. Do not substitute a cumulative mass table containing pre-1980 sedimentation.

## 5. Baseline input to the net-loss calculation

Use `daily_generation.csv` generated by the baseline command and keep its adjacent metadata file. The sediment-hydropower command checks the baseline's parameter/forcing fingerprints and historical period before pairing by `dam_id,date`. Do not compare tables based on their row order, compare different forcing versions, or supply an already sediment-affected baseline.

Output generation is MWh per day; annual generation is the sum of daily MWh. One TWh equals 1,000,000 MWh. Signed net loss is baseline minus sediment-scenario generation.
