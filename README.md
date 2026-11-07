# High Mountain Asia SSY, reservoir sediment and hydropower core models

Core Python research code by **Rui Deng** for specific sediment yield (SSY), reservoir sediment trapping and historical hydropower simulation.

This release contains four workflows:

| Entry point | Purpose |
| --- | --- |
| `scripts/01_train_ssy.py` | Train the 17-predictor, five-member square-root XGBoost model |
| `scripts/02_historical_sediment.py` | Simulate daily reservoir sediment trapping in 1980–2019 |
| `scripts/03_historical_hydropower_baseline.py` | Simulate historical generation with fixed, unsilted active storage |
| `scripts/04_historical_hydropower_sediment.py` | Simulate generation with annual storage updates and calculate signed net differences |

**No confidential observations, reservoir tables, fitted model weights or research results are distributed.** The repository includes source code and code-generated artificial tests only. Plotting, the 12-predictor comparison, calibration and future scenarios are outside this release.

## Scope and reproducibility

The SSY input configuration expects **1,963 records** and 17 predictors. A seeded 70/30 record split produces 1,374 training and 589 validation records. The fixed tree recipe comes from earlier model development; it has not been reoptimized or evaluated on the additional records in this release. No manuscript accuracy scores are encoded or promised.

The reservoir and hydropower workflows start in 1980 and end in 2019 by default. They do not run a pre-1980 spin-up or load an earlier sedimentation state. Original capacities are initialized at the start of the study or commissioning, so the comparison describes **sedimentation accumulated within the simulated historical period**, rather than total sedimentation since construction. It does not reproduce a lifetime calculation that includes an earlier state.

Prepared, locally supplied tables are the starting interface. The source GIS extraction, river-reach matching and environmental-data preparation are not bundled. Users must establish that their prepared inputs use the documented spatial definitions, units and topology.

Read the retained sediment water-volume convention and the annual hydropower timing assumptions in [the methods](docs/methods.md) before interpreting a result. Passing software tests is not independent validation of those scientific assumptions.

## Installation

Use Python **3.11 or 3.12** on a conventional CPU computer. No GPU, cloud account or proprietary software is required. Create an isolated environment, then install from the repository root:

```bash
python -m venv .venv
# Linux/macOS: source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -e ".[test]"
```

Exact core package versions are specified in `requirements.txt` and `pyproject.toml`. Installation time depends on network access and whether wheels are cached. Runtime and memory use for a research run depend on the number of reservoirs and dates; no production benchmark on the confidential data is claimed.

## Check the installation without research data

```bash
python -m pytest
python scripts/smoke_test.py
```

The smoke test creates artificial input tables in a temporary directory, exercises the four workflows, and removes the temporary files. Artificial inputs are not observations, are not fitted to published results, and cannot be used to infer research performance. No dataset file is shipped in the repository.

## Run with locally held inputs

Paths below are placeholders to private files on your own computer. Input schemas, units and missing-value rules are in [input schemas](docs/input_schemas.md).

```bash
python scripts/01_train_ssy.py --input /path/to/private/ssy.csv --output outputs/ssy

python scripts/02_historical_sediment.py --dams /path/to/private/sediment_dams.csv --discharge /path/to/private/discharge.csv --ssy /path/to/private/annual_ssy.csv --output outputs/sediment

python scripts/03_historical_hydropower_baseline.py --dams /path/to/private/hydropower_dams.csv --discharge /path/to/private/discharge.csv --output outputs/baseline

python scripts/04_historical_hydropower_sediment.py --dams /path/to/private/hydropower_dams.csv --discharge /path/to/private/discharge.csv --annual-sediment outputs/sediment/annual_sediment.csv --baseline outputs/baseline/daily_generation.csv --output outputs/sediment_hydropower
```

Use a new output directory for each run. The scripts refuse to overwrite non-empty output directories. Generated model weights and result files are local artifacts and are excluded from version control; do not add them to a public release.

## Documentation

- [Input schemas and units](docs/input_schemas.md)
- [Methods and boundary assumptions](docs/methods.md)
- [SSY core model](docs/ssy.md)
- [Historical sediment kernel](docs/sediment.md)
- [Historical hydropower kernels](docs/hydropower.md)
- [Data availability and confidentiality](DATA_AVAILABILITY.md)
- [Reproducibility and software verification](docs/reproducibility.md)
- [Code availability wording and citation](docs/code_availability.md)

## Citation, versioning and license

The code is released under the [MIT License](LICENSE). Cite a versioned release rather than a moving branch; citation metadata are in [CITATION.cff](CITATION.cff).

Repository: https://github.com/Dr-jf/HMA-SSY-Sediment-Hydropower

Release: https://github.com/Dr-jf/HMA-SSY-Sediment-Hydropower/releases/tag/v1.0.0

A GitHub URL is not a DOI. No DOI or journal endorsement is claimed. For long-term citation, a matching release can also be deposited in a DOI-minting archive such as Zenodo; add the identifier only after the deposit exists.

## Contact

For software questions, open a GitHub issue. Do not post confidential data in issues or pull requests. Questions about access to research data should be directed privately to the corresponding author using the contact information in the associated manuscript.
