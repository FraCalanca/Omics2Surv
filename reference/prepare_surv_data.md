# Prepare Survival Data from MultiAssayExperiment

Extracts survival outcome variables from a MultiAssayExperiment's
colData and optionally reorders samples to match a specified order.

## Usage

``` r
prepare_surv_data(mae, surv_list, sample_order = NULL)
```

## Arguments

- mae:

  A MultiAssayExperiment object containing clinical data in `colData`.

- surv_list:

  Character vector of column names in `colData(mae)` representing
  survival variables (e.g., `c("overall_survival", "status")`).

- sample_order:

  Optional character vector of sample identifiers specifying the desired
  row order. If NULL, original order is preserved. Default is NULL.

## Value

A data frame containing:

- Columns: Survival variables specified in `surv_list`

- Rows: Samples, ordered according to `sample_order` if provided

- Row names: Sample identifiers
