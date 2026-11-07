# Historical hydropower: paired generation and annual capacity updates

The two entry points share one scheduling and energy calculation. The baseline
keeps the initial active-storage capacity fixed. The sediment condition updates
active capacity after generation on the last wet-season day of each year. Their
signed generation difference is computed only after matching every operational
dam/date record.

This is a readable reference implementation of the annual capacity-update
branch. It does **not** implement a later daily-deposition SI formulation. The
distributed checks use artificial values only and do not establish observed
generation accuracy or validate sediment-induced losses.

## Study boundary and fair initialization

The default interval is 1980-01-01 through 2019-12-31. Shorter full-calendar-year
intervals within 1980-2019 are accepted for inspection and tests. No prehistory
or future climate scenarios are simulated.

For each dam, entry is 1 January of `max(start_year, built_year)`. Both conditions
start with the **original active capacity completely filled with water** and
**zero accumulated sediment**. Subsequent days carry actual remaining water;
neither a new year nor a capacity update resets water to full storage. A dam not
yet constructed has no output records, rather than zero-generation observations.

Consequently, a dam built before 1980 has no inherited pre-1980 sediment loss.
These paired runs isolate the effects of deposition accumulated within the
selected study window. They must not be described as reproducing archived
lifetime-accumulation historical indicators. Selecting a later `start_year`
also changes the initialization boundary and resets accumulated sediment there.

## Prepared private inputs

Keep genuine inputs outside the repository. The entry points do not perform
GIS matching, generate discharge, route upstream hydropower releases or infer
missing station parameters.

### Dam parameters (`--dams`)

One unique `dam_id` per row. Identifiers are opaque strings: literal text such
as `NA` is preserved, while empty identifiers are rejected. All four CSV input
roles (parameters, discharge, annual sediment and baseline daily generation)
are checked for duplicate or empty raw header names before parsing; ambiguous
headers are errors rather than silently renamed columns.

| Column | Unit and requirement |
| --- | --- |
| `dam_id` | Nonmissing unique identifier |
| `built_year` | Integer Gregorian construction year; entry is 1 January |
| `active_capacity_m3` | Positive initial active capacity, m³ |
| `dead_capacity_m3` | Nonnegative initial dead capacity, m³ |
| `reservoir_area_km2` | Positive reservoir surface area, km² |
| `dam_height_m` | Positive dam height, m |
| `normal_level_m` | Normal pool elevation, m |
| `dead_level_m` | Dead pool elevation, m |
| `intake_level_m` | Generation intake elevation, m |
| `installed_capacity_mw` | Positive installed-power limit, MW |
| `max_discharge_m3s` | Positive turbine-channel discharge limit, m³/s |
| `wet_months` | Semicolon-separated months, such as `6;7;8;9` |
| `efficiency` | Optional; defaults to 0.90, must be in (0, 1] |

Elevations must use an internally consistent datum and satisfy the estimated
bottom ≤ dead pool < normal pool condition. Wet months must form one continuous
period that does not cross the calendar-year boundary; at least one dry month
must remain. This requirement defines an unambiguous annual update date.

Total initial capacity is active plus dead capacity. With area converted to m²,
the fixed shape parameter and effective intake reference are

\[
m=\operatorname{clip}\!\left(
\frac{A\,H_{\mathrm{dam}}}{C_{A,0}+C_{D,0}},1,10\right),\qquad
z_I^{\ast}=\min(z_I,z_D).
\]

The same fixed values are used in both conditions. The area/height derivation,
bounded exponent and intake-referenced head are idealizations, not measured
storage-elevation geometry or net hydraulic head.

### Daily discharge (`--discharge`)

Columns: `date,dam_id,flow_m3s`. Dates are daily Gregorian dates in 1980-2019;
timezone or subdaily timestamps are rejected. Every selected dam must have
exactly one finite, nonnegative flow for every operational day in the selected
interval. Missing days, duplicate keys and unknown dam IDs are errors. The
forcing may include a longer historical interval than a requested short run;
only operational days in the requested interval are simulated. Preconstruction
values are not used or turned into observations.

### Annual sediment (`--annual-sediment`, sediment condition only)

Columns: `dam_id,year,trapped_t`. `trapped_t` is the **annual incremental trapped
sediment mass in tonnes**, directly compatible with `annual_sediment.csv` from
the historical sediment entry point. Additional columns are ignored.

Every operational dam-year requires one finite nonnegative increment. Annual
values are allowed to rise or fall between years: it is their cumulative sum
that cannot decrease. Other dams and nonoperational years may have missing
values; these are ignored and never interpreted as observed zero deposition.
Duplicate dam/year keys are rejected. Only years in 1980-2019 are accepted.

Do not supply lifetime cumulative mass, Mt, incoming sediment, or a previously
computed capacity-loss percentage under this column name. The adapter cannot
infer the scientific meaning of mislabeled numeric values. Trapping efficiency
and upstream sediment interception must already be represented in `trapped_t`;
they are not applied again by the hydropower model.

## Common daily operation and energy

Let Δt = 86,400 s, `S` be actual active water storage, `C` the capacity applicable
to that day, `Q` inflow and `Qmax` the turbine-channel discharge limit.

During a dry period with `N` days remaining, including today, stored-water
release is

\[
R=\min\left(\frac{S}{N\Delta t},\max(Q_{\max}-Q,0)\right).
\]

Excess inflow is retained if storage space remains:

\[
K=\min\left(\max(Q-Q_{\max},0),\frac{C-S}{\Delta t}\right),\quad
q=\min(Q,Q_{\max})+R,\quad o=\max(Q-Q_{\max}-K,0).
\]

At the start of a wet period the complete specified inflow sequence for that
period is used to plan refill. First use flow exceeding the turbine-channel
limit. If this is insufficient, solve a common additional retention target
subject to each day's available inflow. The target is found by sorting the
available bounded daily flows and distributing the remaining deficit equally
among days that are not already constrained. Only required storage is retained.
If all wet-season inflow is insufficient, refill remains incomplete rather than
creating water. Planning occurs once per wet period in this annual-update
implementation; perfect wet-period inflow knowledge is assumed.

For positive active capacity, water elevation and instantaneous head are

\[
z(S)=z_D+(z_N-z_D)(S/C)^{1/m},\qquad H=z(S)-z_I^{\ast}.
\]

Scheduled discharge is constant within the day. Storage is interpolated linearly
from its start value to its end value before any annual capacity update.
Instantaneous power and daily energy are

\[
P(\tau)=\min\{P_{\max},0.00981\eta q\max[H(\tau),0]\},\qquad
E=\int_0^{24}P(\tau)\,d\tau\quad\mathrm{MWh}.
\]

The power cap is applied **before** integration, including an analytic split
where the changing head crosses the installed-power threshold. The power cap
does not trigger another release optimization. Zero-head routed water is
reported as bypass rather than effective turbine flow.

## Annual sediment allocation and update timing

The bulk density is fixed at 1.2 t/m³. Accumulate only operational annual
increments from study entry:

\[
V_y=\sum_{u=y_{\mathrm{entry}}}^{y}\frac{J_u}{1.2}.
\]

Estimate the bottom as `z_b = z_N - H_dam`, and define

\[
f_D=\frac{z_D-z_b}{z_N-z_b},\qquad f_A=1-f_D.
\]

For each year's cumulative `V_y`, simultaneously allocate `f_D V_y` to initial
dead capacity and `f_A V_y` to initial active capacity, each capped by its zone's
capacity. Reassign excess first to available active capacity, then to available
dead capacity. Log material exceeding both zones as unallocated volume. Recompute
this allocation from the cumulative volume once; do not add cumulative totals
repeatedly. The remaining active capacity is

\[
C_{A,y}=\max(0,C_{A,0}-D_{A,y}).
\]

The annual update happens **after generation on the final wet-month day**.
For example, a June-September wet period updates after 30 September generation,
with the reduced capacity used on 1 October. The **entire calendar year's**
sediment is concentrated at this node, including sediment assigned to later
months. This is an explicit annual approximation, not measured deposition up
to that wet-end date and not daily deposition at the beginning of each day.

Normal/dead levels, intake, `m`, efficiency, installed power and discharge limit
remain fixed. Capacity changes compress the horizontal storage scale of the
fixed elevation relation. They do not reconstruct a sediment-modified E-V-A
curve. Only active water storage is carried by this branch; dead-capacity loss
is tracked as sediment allocation and does not create a separate transfer of
dead-pool water into active water.

At an update, clip actual water to the new capacity:

\[
S^+=\min(S^-,C_{A,y}),\qquad O_{\mathrm{capacity}}=S^- - S^+.
\]

This excess water is **non-generating spill**, included in the day's water
balance. For all days the active-water ledger satisfies

\[
S_{\mathrm{end}}-S_{\mathrm{start}}
=\Delta t(Q-q)-O_{\mathrm{hydraulic}}-O_{\mathrm{capacity}}.
\]

The calculation checks storage/capacity bounds, nonnegative flows, discharge
and power limits, sediment allocation and active-water closure. Numerical
water tolerance is `max(1e-5, 1e-12 * relevant_volume_scale)` m³; it does not
constitute physical validation.

When active capacity reaches exactly zero, set active water to zero and stop
refill/release. Prescribe normal-pool head `z_N - z_I*`, route `min(Q,Qmax)` and
spill the remaining inflow. Do not divide storage by zero. This idealized
run-of-river boundary assumes persistent hydraulic access and does not predict
plant survival after sedimentation.

## Signed loss, output and provenance

The sediment entry point calculates

\[
L_t=E_t^0-E_t^s,\qquad L_y=\sum_{t\in y}L_t.
\]

Negative daily differences are retained: changed storage and refill timing can
increase generation on some days. No seasonal zero-loss rule is imposed. Annual
loss fractions use summed net loss divided by summed baseline generation;
a zero baseline denominator is undefined. Cumulative sediment/capacity states
must not be summed over days.

Each output directory contains:

- `daily_generation.csv`: operational dam/date rows, daily energy (MWh), flow
  (m³/s), water/capacity/volume-ledger fields (m³), mean head (m), and numerical
  residuals. The sediment output additionally has `baseline_energy_mwh` and
  signed `net_loss_mwh`; its `energy_mwh` is sediment-condition generation.
- `annual_generation.csv`: calendar-year energy and spill sums; sediment runs
  additionally report annual baseline energy, signed net loss and loss fraction.
- `manifest.json`: completion status, method and implementation fingerprints,
  interval, input/output SHA256 hashes and row count. No local input paths are
  recorded.

`routed_m3s` is all scheduled hydraulic outflow through the turbine/bypass route;
`turbine_m3s + bypass_m3s` equals it. `spill_m3` includes hydraulic spill and
capacity-update spill. On update days, energy/head and start capacity refer to
operation before the update, while end capacity/storage refer to its result.

For paired CLI runs the baseline daily CSV requires its adjacent `manifest.json`.
The sediment command checks identical dam/discharge content hashes, study
interval, implementation fingerprint and method version, plus baseline-file
integrity. It then rejects differing or duplicate daily keys rather than taking
an inner join that discards records. The in-memory `calculate_net_loss` function
checks keys and energy; its caller remains responsible for identical forcing,
parameters and initialization. SHA256 records provide consistency checks, not
authenticated claims about externally supplied data.

Output directories must be new or empty. Existing results are never overwritten.
Use a directory only after its manifest reports `COMPLETED`.

## Commands and Python API

From the repository root, with privately prepared inputs:

```text
python scripts/03_historical_hydropower_baseline.py --dams private/hp_dams.csv --discharge private/discharge.csv --output outputs/baseline
python scripts/04_historical_hydropower_sediment.py --dams private/hp_dams.csv --discharge private/discharge.csv --annual-sediment outputs/sediment/annual_sediment.csv --baseline outputs/baseline/daily_generation.csv --output outputs/hydropower_sediment
```

Both accept `--start-year` and `--end-year` within 1980-2019, inclusive. The
defaults are 1980 and 2019. The two runs must use the same interval and exact
same dam/discharge files.

The core exports `HydropowerDam`, `load_dams_csv`, `simulate_baseline`,
`simulate_sediment` and `calculate_net_loss`. Both simulators return daily
pandas DataFrames. `daily_energy_mwh` and `allocate_cumulative_sediment` expose
the energy and allocation equations for hand-checkable scientific inspection.

This model excludes turbine abrasion, sediment-removal adaptation, evaporation,
seepage, non-hydropower demands, grid constraints, outages, pumping consumption,
tailwater variability and coordinated cascade operation. No observed-loss
validation, uncertainty propagation or reproduction of confidential study
outputs is claimed.
