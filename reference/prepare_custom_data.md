# Prepare Custom Data in LinkedOmics-like Format

Validates and cleans a user-supplied omics or clinical matrix, applying
the same quality checks and NA handling used by
[`retrieve_linkedomics_data`](https://fracalanca.github.io/Omics2Surv/reference/retrieve_linkedomics_data.md).
Use this function when data come from a source other than LinkedOmics
but must conform to the same format expected by the rest of the
pipeline.

## Usage

``` r
prepare_custom_data(data, data_type = c("Omics", "Clinical"))
```

## Arguments

- data:

  A data frame or matrix in LinkedOmics-like format: features as rows,
  samples as columns, with row names set to feature identifiers.

- data_type:

  Character string specifying the data type:

  - `"Omics"`: features with any `NA` are removed.

  - `"Clinical"`: samples with `NA` in `"overall_survival"` or
    `"status"` are removed.

## Value

A matrix with:

- For `"Omics"`: features (rows) with any `NA` removed.

- For `"Clinical"`: samples (columns) with `NA` in survival variables
  removed.

## See also

[`retrieve_linkedomics_data`](https://fracalanca.github.io/Omics2Surv/reference/retrieve_linkedomics_data.md)
