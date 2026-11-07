# Historical reservoir sediment model

Author: Rui Deng

This readable reference implementation accepts prepared tables and models a
complete-calendar-year interval within **1980–2019**. The default is the full
interval; a shorter interval supports artificial examples. There is no
pre-study simulation, future scenario, spatial extraction, sediment management,
calibration, or bundled research dataset. The study's inventory size is not an
API restriction. Tests use invented identifiers and quantities generated in code.

## Prepared inputs

All identifiers are opaque, nonempty strings. A single, consistent inventory is
required across the tables; input row order is irrelevant. Raw CSV column
headers must be nonempty and unique. The CLI validates them before pandas can
rename duplicate fields; literal identifiers such as `NA` and leading-zero
strings are preserved rather than converted to missing values or numbers.

| Table | Required columns | Units or definition |
|---|---|---|
| Dams | `dam_id`, `built_year`, `total_capacity_m3`, `dead_capacity_m3`, `basin_area_km2`, `upstream_ids` | Calendar year; capacities in m³; **full upstream control-basin** area in km²; direct upstream IDs separated by semicolons in CSV, empty for a headwater dam |
| Discharge | `date`, `dam_id`, `flow_m3s` | Exactly one timezone-free calendar-day row for every dam in the selected interval, including leap days; m³/s |
| Annual SSY | `year`, `dam_id`, `ssy_t_km2_yr` | Exactly one row per dam and selected year; full-control-basin mean t km⁻² yr⁻¹ |

The program rejects missing/duplicate keys, incomplete calendars, extraneous
dates/years/IDs, fractional years, negative/nonfinite forcing, nonpositive total
capacity or basin area, dead storage outside total capacity, unknown or repeated
upstream IDs, and network cycles. Each dam may feed at most one downstream dam:
no unparameterized bifurcation may duplicate its exported sediment. Prepared
data are required even before commissioning because the unregulated sediment
still passes downstream. No missing parameter is silently estimated.

`upstream_ids` must identify **immediate** upstream model nodes, not a list of all
ancestors. The prepared interface trusts the supplied basin boundaries, areas,
and network. It does not reconstruct or verify their geographic validity.

## Initialization and interpretation

Every reservoir starts with its original total capacity. On its first operating
day, stored water is initialized to reported dead capacity and suspended mass
to zero. For a dam completed before the selected start year, that first day is
1 January of the start year. For a dam completed during the interval, it is
1 January of its commissioning year. The initial water is an explicit boundary
addition, recorded as `initial_water_added_m3`.

The 30-day storage window is seeded with that initial dead-pool water; its
outflow window is seeded with mean direct discharge in the **first selected
study year** × 86400. This same first-year boundary is used for later-commissioned
dams; no unobserved earlier state or earlier year is read. These are idealized
initial conditions, not reconstructed historical reservoir states.

Consequently, cumulative sediment quantities cover only the selected study
period after commissioning. They are **not lifetime sedimentation**, and this
release cannot reproduce a run that inherits pre-study deposition. There is no
annual reset: water, suspended mass, capacity and the rolling windows carry
across years. A hypothetical known start-of-study sediment state is not invented
or inferred by this API.

Before commissioning, water and sediment pass through. `is_operating` is false,
and trapping, remaining-capacity, storage-state and TE fields are missing rather
than presented as zero-valued observations. Both operating and inactive rows
are model records, never observations.

## Annual input and daily network calculation

For full-basin SSY \(Y_{i,y}\), basin area \(A_i\), and immediate upstream set
\(U_i\), annual masses are

\[
G_{i,y}=Y_{i,y}A_i,\qquad
L_{i,y}=\max\left(G_{i,y}-\sum_{j\in U_i}G_{j,y},0\right).
\]

`negative_increment_clipped` exposes any negative gross difference. Clipping is
a modeling choice and can indicate inconsistent spatial forcing; it is not an
independent validation. Local annual mass is divided among the year's days with
weights \(Q_{i,d}^2/\sum_d Q_{i,d}^2\). Scaling Q by its annual maximum before
squaring preserves those weights and avoids numerical overflow. If all Q is zero
and local annual mass is positive, the retained legacy rule divides mass
uniformly over the year and sets `uniform_mass_fallback`. This fallback is not a
physically justified dry sediment transport model.

Every day is calculated in upstream-to-downstream topological order. Incoming
sediment is local daily mass plus the **same day's** direct upstream exported
mass. Water is independently forced at each dam:
\(I_{i,d}=86400\max(Q_{i,d},0)\). There is no basin-area scaling of Q, upstream
water addition, water routing, or channel-travel-time lag.

For the preceding 30 stored-water and actual daily-outflow values, the daily
residence time and Brune/Gill middle-curve efficiency are

\[
\tau=\frac{\overline S}{365.25\overline O},\qquad
TE=\frac{\tau}{0.012+1.02\tau}.
\]

The current day's results are inserted into the windows **after** calculating
TE. If mean outflow is zero, tau is \(10^9\) years for positive mean storage and
zero otherwise. TE is limited to [0,1], and is zero once capacity is exhausted.
Dry bulk density is fixed at \(\rho=1.2\) t/m³. For incoming mass \(M\), old
suspended mass \(B\), remaining capacity \(C\), and old stored water \(S\):

\[
D=\min(TE\,M,\rho C),\quad v=D/\rho,\quad C'=\max(C-v,0).
\]

Only **today's incoming mass**, not the old suspended stock, enters the
deposition formula. The old stock enters the subsequent mixed export balance.
Capacity limits deposition; there is no dredging, flushing, resuspension,
compaction, evaporation, operational release target or reservoir-specific
sediment-density calibration.

## Explicit retained water-volume convention

The implementation deliberately retains the archived research equation:

\[
W=\max(S+I-v,0),\quad O=\max(W-C',0),\quad
S'=\min(C',\max(W-O,0)).
\]

Thus deposited bulk volume both reduces capacity and is subtracted from water.
It is an explicit **water-volume sink convention**, not a newly corrected
physical displacement calculation. This code does not claim physical water
conservation, site-specific physical validation, or a separate pore-water model.
Revising that convention would change the algorithm and its resulting TE and
would require a separately documented scientific decision.

Remaining sediment \(R=\max(B+M-D,0)\) is mixed using concentration
\(c=R/(S+I)\), or zero for a zero denominator. Export is
\(M_{out}=\min(R,cO)\); the rest is the next suspended stock. Each day reports

\[
r_w=S_{opening}+S_{initial}+I-v-O-S',\qquad
r_s=B+M-D-M_{out}-B'.
\]

`water_sink_m3` equals v. The floor at zero in W can leave a nonzero water
residual if the requested sink exceeds available water; this is retained and
exposed, not forced to close. A small residual only checks the specified
internal arithmetic. It does not validate the research convention or the input
data. The artificial tests explicitly include this conflict.

## API and command line

```python
from hma_models.sediment import (
    SedimentDam, simulate_historical_sediment, summarize_annual,
)

# Private prepared DataFrames and parameters are supplied by the caller.
daily = simulate_historical_sediment(
    dams, discharge, annual_ssy, start_year=1980, end_year=2019,
)
annual = summarize_annual(daily)
```

`dams` is a sequence of `SedimentDam` records; Python `upstream_ids` is a tuple
or list of strings. The other two inputs are pandas DataFrames following the
schemas above.

```console
python scripts/02_historical_sediment.py --dams private/dams.csv --discharge private/discharge.csv --ssy private/annual_ssy.csv --output outputs/sediment --start-year 1980 --end-year 2019
```

The destination must be new or empty. The CLI writes `daily_sediment.csv` and
`annual_sediment.csv` locally and refuses a nonempty destination. Input and
generated-output files are not included in this repository or uploaded by the
program.

## Output fields

Daily rows have `date`, `year`, `dam_id`, `is_operating`, and:

| Fields | Definition |
|---|---|
| `local_sediment_t`, `incoming_sediment_t`, `outgoing_sediment_t` | Local input and total network input/export, tonnes per day |
| `trapped_t`, `storage_loss_m3` | New daily deposited mass and volume, t/day and m³/day |
| `cumulative_trapped_t` | Cumulative mass since the selected study/commissioning start, tonnes |
| `remaining_capacity_m3`, `stored_water_m3`, `suspended_t` | End-of-day capacity, water and suspended state |
| `trapping_efficiency`, `residence_time_years` | Prior-window middle-curve TE, dimensionless; tau in years |
| `water_in_m3`, `water_out_m3`, `initial_water_added_m3`, `water_sink_m3` | Daily water ledger, m³ |
| `water_balance_residual_m3`, `sediment_balance_residual_t` | Daily specified-convention ledger residuals, m³ and tonnes |
| `negative_increment_clipped`, `uniform_mass_fallback` | Explicit annual-forcing flags repeated on affected dam-days |

Annual columns are `year`, `dam_id`, `trapped_t`, `storage_loss_m3`,
`cumulative_trapped_t`, `is_operating`. Annual amounts sum daily increments;
cumulative mass is the year-end study-period total. An entirely uncommissioned
year remains missing for these amounts. Downstream hydropower code consumes
**annual new `trapped_t`**, not the cumulative column, to avoid double counting.
An annual deposited-mass/incoming-mass trapping ratio is not supplied or
interchanged with TE. These outputs do not establish independent validation
against observed annual sediment retention.
