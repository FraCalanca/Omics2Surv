# Prepare Design Matrix from MultiAssayExperiment for Modeling

Extracts omics data from a MultiAssayExperiment object and optionally
combines it with clinical variables to create a design matrix suitable
for survival modeling. Handles sample ordering, NA removal, and factor
variable conversion.

## Usage

``` r
prepare_design_matrix(mae, assay_name, clinical_list = NULL, remove_na = TRUE)
```

## Arguments

- mae:

  A MultiAssayExperiment object containing omics experiments and
  clinical data in `colData`.

- assay_name:

  Character string specifying the name of the assay (omics layer) to
  extract from the MAE. Must match an experiment name in `names(mae)`.

- clinical_list:

  Character vector of clinical variable names to include in the design
  matrix. Variables are extracted from `colData(mae)`. If NULL, only
  omics features are included. Default is NULL.

- remove_na:

  Logical indicating whether to remove samples (rows) with any missing
  values. If TRUE, samples with NA in any variable are excluded. If
  FALSE, NAs are retained. Default is TRUE.

## Value

A numeric matrix where:

- Rows represent samples (patients)

- Columns represent features:

  - Clinical variables (if `clinical_list` is specified), followed by

  - Omics features (genes, proteins, metabolites, etc.)

- Row names are sample identifiers

- Column names are feature names

- All values are numeric (factors are converted to numeric codes)

## Details

The function performs the following operations:

1.  **Extract omics data**: Retrieves the "counts" assay from the
    specified experiment in the MAE

2.  **Determine sample order**: Uses column names from omics data as the
    reference sample ordering

3.  **Extract clinical data**: Retrieves clinical variables from MAE
    colData, ordered to match omics samples

4.  **Combine data** (if clinical_list is specified):

    - Subsets clinical variables to requested list

    - Converts factors to numeric codes

    - Binds clinical variables (columns) with transposed omics data

5.  **Handle missing data** (if remove_na = TRUE):

    - Identifies samples with any NA values

    - Removes these samples from the design matrix

## See also

[`MultiAssayExperiment`](https://github.com/waldronlab/MultiAssayExperiment/reference/MultiAssayExperiment.html)
for the MAE object class `colData` for accessing clinical data
[`prepare_surv_data`](https://fracalanca.github.io/Omics2Surv/reference/prepare_surv_data.md)
for preparing survival outcome data
