# Data availability and confidentiality

Research observations and prepared inputs are confidential and are **not included** in this repository. This applies to SSY observations, environmental predictor tables, station identifiers, reservoir attributes, discharge series, sediment tables and fitted model weights.

The default SSY input count is 1,963 records, as specified for this release. No private records have been appended, reconstructed, published or used to rerun the research during preparation of the release.

The code accepts user-supplied local tables with the schemas in `docs/input_schemas.md`. The schemas describe column names and units, not measured values. Unit tests and the smoke test create explicitly artificial numbers in memory or temporary files; no observational or surrogate research dataset is distributed.

## Access to research data

Data access is governed by the relevant data owners and confidentiality restrictions. Requests should be made privately to the corresponding author through the contact information in the associated manuscript. This code release does not grant access, promise that every request can be approved, or establish terms or response times on behalf of a data owner.

For journal submission, the manuscript's Data availability statement must identify the actual reason for each restriction and the applicable access procedure. Those details, including any terms, contacts and editorial arrangements, must be supplied by the authors and data owners; they cannot be inferred from a public code repository. Do not claim that a code-only release alone establishes compliance with data-sharing requirements.

## Keeping generated artifacts private

Keep research inputs outside the checkout. Write outputs under a local ignored directory or another private location. The `.gitignore` excludes common data and model-output locations, but an ignore file does not sanitize files already committed. Inspect the staged file list before every upload. Never include confidential data, fitted tree files, predictions, notebook outputs or credentials in a release or issue attachment.
