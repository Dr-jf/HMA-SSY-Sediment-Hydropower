# Methods and boundaries

## SSY

The core model uses 17 ordered environmental predictors and a square-root target. Five XGBoost members use the fixed recipe documented in `ssy.md`. The first member uses all training records and the other four use seeded, non-replacement 85% samples. Predictions are averaged on the square-root scale, clipped below at zero, then squared. No affine calibration, regional offset, plotting or uncertainty interval is added.

The seeded 70/30 split is by records, not by station or forward time. It does not establish independent spatial or future-climate generalization. The configured 1,963-record count is not evidence of fitting or validation on those private records.

## Historical sediment

The kernel uses supplied upstream relationships, direct independent discharge at each reservoir, incremental annual sediment yield and Q-squared daily allocation. Only sediment is cascaded between reservoirs. A reservoir's water forcing is not increased by upstream reservoir outflow. Daily trapping uses the preceding 30-day storage/outflow window and the Brune-Gill M curve, with density 1.2 tonnes per cubic metre.

This reference preserves the archived **sediment-volume deduction from stored water** convention: deposited volume reduces both remaining capacity and the water quantity used before release. Its internal water ledger therefore includes a deposited-volume sink. This is an explicit research convention retained for algorithm transparency, not a claim that the conventional physical water balance has been independently validated or repaired. Changing it changes subsequent outflow, residence time and trapping; such a change requires a separately versioned method and new scientific runs.

Old suspended sediment participates in mixing and export but is not added to the same day's trapping formula. The kernel does not include erosion, flushing, re-suspension, compaction or operational sediment management. Supplied topology must distinguish complete upstream catchments from local incremental sediment sources.

## Historical initialization

There is no pre-1980 calculation. For a previously commissioned reservoir, the start of the study uses original capacity, dead-pool water and zero suspended sediment in the sediment model. A reservoir commissioned during the study activates on 1 January of its commissioning year. The initial 30-day outflow buffer uses the mean discharge of the first study year. These are explicit initialization assumptions, not recovered observations of the 1980 state.

Both hydropower scenarios instead initialize the **active-storage water state** at full capacity on the same study/commissioning date. Sediment cumulative volume starts at zero. Neither model silently imports a lifetime cumulative sediment table.

## Annual hydropower storage update

Baseline active capacity is fixed. In the sediment scenario, annual incremental trapped mass is accumulated from the historical simulation start, divided by 1.2 t/m3, and allocated between dead and active storage using the prescribed height fractions, with overflow reassigned to available space.

On the last day of the reservoir's final wet month, generation is first calculated with the old capacity. The active capacity is then updated using the **whole calendar year's** cumulative sediment. The next day uses the new capacity. Thus sediment deposited after the wet-season endpoint is represented earlier than its physical deposition date; this is the selected annual approximation. Stored water above the new capacity becomes non-generating spill and is included in the water balance.

Only active capacity changes in the operating geometry. Shape exponent, pool levels, intake elevation, efficiency, installed capacity and turbine discharge limit remain fixed. When active capacity reaches zero, the stated idealized boundary is run-of-river generation at normal-pool head. This is not a daily sediment-geometry reconstruction.

Daily net generation change is baseline generation minus sediment-scenario generation. Negative differences are retained and enter annual net totals. Counterfactuals must share the same dams, forcing, dates and fixed parameters; the command-line workflow checks their input fingerprints.
