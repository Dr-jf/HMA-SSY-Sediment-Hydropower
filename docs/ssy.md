# SSY machine-learning core

Author: Rui Deng.

This module implements the previously selected 17-predictor, five-member
square-root XGBoost recipe through a private tabular-input interface. The default
1,963 records are user-reported (1,911 plus 52); this repository does not contain
or independently verify those observations. The inherited parameters have not
been retuned or validated on that enlarged sample. Artificial tests establish
software behavior, not manuscript performance.

## Private input schema

Supply a UTF-8 CSV with one station-year record per row, a header, the 17 columns
below and the target `SSY`. The program selects predictors in this fixed order
even if CSV columns are reordered. The delta-symbol aliases `ΔG` and `ΔSC` are
accepted in place of `DeltaG` and `DeltaSC`; supplying both names for the same
predictor is rejected. Extra columns are ignored and are never saved.

| Column | Meaning and required input units | Monotonic direction |
|---|---|---:|
| PGA | Peak ground acceleration in the original source's numeric unit | +1 |
| S | Mean slope, degrees | +1 |
| V | NDVI, dimensionless | -1 |
| P | Annual precipitation, mm/year | +1 |
| T | Annual mean temperature, degrees Celsius | +1 |
| L | Lithologic erodibility, dimensionless source index | -1 |
| Ksn | Channel steepness index in the original source's numeric unit | +1 |
| A | Drainage area used as the predictor, km² | -1 |
| G | Glacier cover, percent | +1 |
| F | Permafrost cover, percent | +1 |
| R | Elevation relief, maximum minus minimum, m | +1 |
| GMC | Glacier mass change, m water equivalent/year | -1 |
| Q | Annual surface runoff, mm/year; not discharge in m³/s | +1 |
| SC | Snow cover, percent | +1 |
| SM | Annual snowmelt total, mm water equivalent/year | +1 |
| DeltaG | Annual glacier-cover change, percentage points/year | -1 |
| DeltaSC | Annual snow-cover change, percentage points/year | -1 |
| SSY | Observed annual specific sediment yield, t km⁻² year⁻¹ | target |

PGA and Ksn retain source-dependent numeric units; no conversion is invented
here. Match the upstream feature definitions and units before training. This
module neither derives environmental predictors nor harmonizes their sources.
Predictors may contain missing values (`NaN` or empty CSV cells), which follow
XGBoost's native missing-value branches. Infinite or nonnumeric predictor values
are rejected. `SSY` must be finite and nonnegative, with no missing targets.
Duplicate column names and missing required columns are rejected. There is no
imputation, feature scaling or target adjustment.

## Training sequence

1. Check the schema, numeric values and exact expected number of rows.
2. Use a shuffled 70:30 row split with seed 42. For 1,963 rows, the split contains
   1,374 training and 589 held-out rows. Counts are computed from the supplied
   sample; the old sample's partition sizes are not hardcoded. Preserve the
   input row order to reproduce this positional split.
3. Fit five members to the square root of SSY. Member one uses all training
   rows; members two through five each use a deterministic, without-replacement
   85% sample. Their subset seeds are `20261008 + 5000 + j`, for zero-based member
   indices `j = 1, 2, 3, 4`. Sample size is `max(2, floor(0.85 * n_training))`.
4. Use fixed tree counts **788, 736, 465, 316, 330**, the same 17 monotonic
   directions, `reg:squarederror`, `hist`, and model seed 20261008. No held-out
   responses are passed to fitting. There is no early stopping or parameter
   search in this core: tree counts and parameters were selected during earlier
   model development.
5. Save five native JSON models and aggregate metadata. `predict_ssy` loads these
   weights and returns:

   `SSY_raw = max(mean(member_1, member_2, member_3, member_4, member_5), 0)²`

The mean is taken on native square-root predictions, before one nonnegative
square. There is no affine calibration, scenario offset, uncertainty analysis,
predictor-removal branch or plotting step. The held-out partition is reserved;
this core deliberately does not export per-record observations or predictions.

The inherited base parameters are:

| Parameter | Value |
|---|---:|
| max_depth | 6 |
| learning_rate | 0.02262021088229677 |
| min_child_weight | 0.06737187721626481 |
| subsample | 1.0 |
| colsample_bytree | 0.9779304618759297 |
| reg_lambda | 0.01032645739707662 |
| reg_alpha | 1.0 |
| gamma | 0.0019485455718767468 |
| max_bin | 512 |
| grow_policy | depthwise |

The member-level 85% sample is separate from XGBoost's `subsample=1.0`.

## CLI and Python API

After installing the project dependencies, run from the repository root:

```sh
python scripts/01_train_ssy.py --input private_inputs/ssy.csv --output outputs/ssy_model
```

`--expected-samples` defaults to 1963. An explicit alternative permits another
declared input size; it does not verify the provenance of those rows.
`--split-seed` defaults to 42 and `--threads` defaults to 1. All five members are
fitted sequentially. A nonempty output directory is rejected, preserving its
contents. A failed training run can leave partial files; use a fresh destination
for another attempt. Do not publish `outputs/` or any fitted weights.

```python
import pandas as pd
from hma_models.ssy import SSYConfig, predict_ssy, train_ssy

frame = pd.read_csv("private_inputs/ssy.csv", float_precision="round_trip")
metadata = train_ssy(frame, "outputs/ssy_model", config=SSYConfig())
predictions = predict_ssy("outputs/ssy_model", frame)
```

The prediction interface requires the 17 predictors, but not `SSY`; it accepts
any number of prediction rows. It checks model hashes and feature order before
inference. Loading artifacts from an unrelated model is not supported. Neither
API returns a validation-performance claim.

`SSYConfig(expected_samples=20, member_trees=(2, 2, 2, 2, 2))` is intended for a
small, code-generated artificial smoke test. Such changed tree counts are not
the full inherited recipe. `model_seed` is also explicit in the Python API;
changing seeds changes the experiment.

The output consists only of `member_01.json` through `member_05.json` and
`metadata.json` (parameters, versions, aggregate counts, seeds and hashes).
No input paths, station identifiers, row identifiers, original observations or
per-record predictions are saved. Fitted weights can nevertheless encode
private training information and must remain local and untracked.

## Interpretation limits

The split is random by row, not by station or time. It does not establish
independent spatial, temporal or external validation. Previous selection of the
recipe is not repeated here. Results depend on private inputs, their row order,
numeric precision and software versions. This pure core does not reproduce
figures that used a later calibrated model, and no previous manuscript metric
is embedded as an expected result. Reproduction of confidential numerical
results requires the authorized source data and the appropriate complete
research workflow, which are not distributed in this repository.
