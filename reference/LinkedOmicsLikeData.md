# LinkedOmics-like data format

This package accepts either data retrieved from LinkedOmics or
user-provided data, as long as the input follows the same structure as
LinkedOmics data after import.

## Details

The expected input is a feature-by-sample matrix.

**General requirements:**

- Input must be a `data.frame` or `matrix`.

- Rows represent features (genes, proteins, CpGs, or clinical
  variables).

- Columns represent samples.

- Feature identifiers must be stored as row names.

- Sample identifiers must be stored as column names.

**Omics data:**

- Each row corresponds to a molecular feature.

- Values should be numeric.

**Clinical data:**

- Rows correspond to clinical variables.

- The following rows are required:

  - `overall_survival`

  - `status`
