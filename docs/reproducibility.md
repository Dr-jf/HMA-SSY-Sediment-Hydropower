# Reproducibility and software verification

## What is reproducible from this release

This is a source-code release with explicit private-input interfaces. A person with suitable authorized inputs can run the four documented core workflows. Someone without the confidential inputs can install the package, inspect the algorithms, execute artificial unit tests and run the artificial end-to-end smoke check.

The release does not include observations, prepared basin data, trained research trees, predictions or production results. It does not claim to reproduce manuscript numbers without those inputs. It omits plotting, affine calibration, the reduced-predictor comparison, future scenarios and pre-1980 state generation, as stated in the README.

## Environment and installation

Python 3.11 and 3.12 are the supported interpreter range. The core packages are pinned in `requirements.txt` and `pyproject.toml`. `requirements-lock.txt` records the complete package set observed in the tested Windows Python 3.12 environment; other systems can need additional platform-specific dependencies.

Standard installation:

```bash
python -m venv .venv
# Activate the environment for your operating system.
python -m pip install -e ".[test]"
python -m pip check
```

For the recorded dependency set:

```bash
python -m pip install -r requirements-lock.txt
python -m pip install --no-deps -e .
```

No non-standard hardware is required for the artificial checks. Installation requires network access unless the dependencies are cached; allow several minutes as a planning estimate, not a measured universal installation time. Full-run time and memory depend on input size. No confidential-data runtime benchmark is supplied.

## Tests and artificial demo

```bash
python -m pytest
python scripts/smoke_test.py
```

Release preparation on Windows with Python 3.12 used a newly created environment and the pinned core dependencies. **86 unit and command-line tests passed.** Observed full-suite elapsed times were on the order of tens of seconds. The four-entry-point smoke check completed in approximately **6–7 seconds** on the preparation machine. These are illustrative observed timings, not performance promises.

The smoke check creates 40 artificial SSY records, a single artificial reservoir and one leap-year discharge sequence in a temporary directory. It runs each real command-line entry point, passes annual incremental sediment directly into the sediment-hydropower command, checks its pairing against the baseline and removes the temporary inputs and outputs. No artificial dataset is committed as a file, and no research-data computation is performed.

The tests exercise:

- the 1,963-record count and dynamic 1,374/589 row split;
- held-out target isolation, feature order, delta aliases, missing predictor routing and model-file integrity;
- CSV header ambiguity, string identifiers, duplicate/missing dates and invalid values;
- topology, sediment mass accounting, previous-only trapping windows and commissioning boundaries;
- fixed and annually updated water storage, capacity-update spill, power bounds and signed paired differences;
- rejection of incompatible baseline provenance and non-empty output destinations.

Artificial checks establish these software properties. They do not independently validate the model's environmental relationships, reservoir operating assumptions or sediment-volume water convention. Refer to `methods.md` for those limits.

## Determinism and version records

SSY splitting and member subsampling use fixed seeds. The model's metadata records the predictor order, package versions, tree recipe, row counts and member hashes. Keep the input row order fixed if reproducing a positional split. Numerical library versions and platforms can still affect exact floating-point outputs.

Hydropower manifests bind the period, parameter file, discharge file, implementation and saved baseline output using SHA-256. The sediment-adjusted command checks these fingerprints before comparing generation. Use new output directories for changed inputs; do not bypass a mismatch or relabel one method as another.

## Publication archive and privacy

Use the versioned GitHub release link and `CITATION.cff` to identify the source version. `SHA256SUMS` lists the published text-file fingerprints. Generated data and weights remain outside the public source archive. A DOI is available only if a separate DOI-minting archive actually records the release; no such identifier is asserted here.

GitHub Actions defines the same artificial checks for Python 3.11 and 3.12 on Linux. Its actual run status is visible in the repository's Actions tab; a workflow file by itself is not evidence of a successful remote run.
