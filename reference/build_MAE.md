# Build a MultiAssayExperiment Object from Omics Data

Constructs a MultiAssayExperiment (MAE) container that integrates
multiple omics datasets with associated clinical/phenotypic data. This
unified data structure facilitates multi-omics analysis and ensures
coordinated sample management across different data types.

## Usage

``` r
build_MAE(omics_list, clinical_data)
```

## Arguments

- omics_list:

  A named list of SummarizedExperiment objects, where each element
  represents a different omics dataset (e.g., gene expression,
  methylation, proteomics). All elements must be valid
  SummarizedExperiment objects. List names will be used as experiment
  identifiers in the resulting MAE.

- clinical_data:

  A matrix or data frame of clinical/phenotypic data where:

  - Rows represent different clinical variables (features)

  - Columns represent samples

  - Column names should match sample identifiers across omics datasets

  This will be transposed internally so that rows represent samples in
  the final MAE.

## Value

A MultiAssayExperiment object containing:

- `experiments`: The omics datasets from `omics_list`

- `colData`: Transposed and type-converted clinical data as a DataFrame,
  where rows represent samples and columns represent clinical variables

- Coordinated sample mapping across all experiments

## Details

The function performs the following operations:

1.  Validates that all elements in `omics_list` are SummarizedExperiment
    objects

2.  Transposes `clinical_data` so samples are in rows (standard MAE
    format)

3.  Preserves sample identifiers from rownames after transposition

4.  Automatically converts column types using `type.convert` (e.g.,
    strings to factors, numeric strings to numbers) with `as.is = FALSE`

5.  Constructs the MultiAssayExperiment with proper S4 DataFrame
    structure

6.  Prints a confirmation message with the number of omics datasets
    included

The resulting MAE object enables:

- Unified access to multiple omics layers

- Automatic subsetting across all experiments

- Integration with Bioconductor workflows

- Sample-level metadata management

## Data Structure Requirements

- Each SummarizedExperiment in `omics_list` should have sample
  identifiers as column names (in `colData`)

- Clinical data column names should correspond to these sample
  identifiers

- Sample matching across experiments is handled automatically by MAE

## See also

[`MultiAssayExperiment`](https://github.com/waldronlab/MultiAssayExperiment/reference/MultiAssayExperiment.html)
for the MAE class structure
[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
for omics data structure
