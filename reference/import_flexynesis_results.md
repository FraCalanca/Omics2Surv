# Import flexynesis predictions

Reads a `predicted_labels.csv` and returns one prediction per sample for
a single split. The file holds one row per (sample, variable) for
numeric targets and per (sample, variable, class) for categorical
targets, so a naive `row.names = sample_id` (as in the original code)
would hit duplicate row names. This version filters to a single
`variable`, de-duplicates by sample, and (optionally) verifies that no
expected patient is missing.

## Usage

``` r
import_flexynesis_results(
  csv_path,
  split = c("test", "train"),
  variable = NULL,
  expected_samples = NULL,
  numeric_pred = TRUE
)
```

## Arguments

- csv_path:

  Path to a flexynesis `predicted_labels.csv`.

- split:

  Which split to extract: `"test"` (default) or `"train"`.

- variable:

  Optional outcome/variable name to filter on. If `NULL` and several are
  present, the first is used (with a warning).

- expected_samples:

  Optional character vector of sample IDs that must all be present;
  missing IDs raise an error (patient-loss guard).

- numeric_pred:

  Logical; coerce predictions to numeric (TRUE for survival risk scores
  / regression). Default `TRUE`.

## Value

A one-column data frame of predictions (`pred`) indexed by sample
identifier, ordered to match `expected_samples` when supplied.
