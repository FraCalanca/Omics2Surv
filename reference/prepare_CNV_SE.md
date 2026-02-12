# Prepare Copy Number Variation Data as SummarizedExperiment

Processes copy number variation (CNV) data by filtering constant genes,
optionally performing feature selection using supervised
(survival-based) or unsupervised (variability-based) methods, and
packaging the results into a SummarizedExperiment object for downstream
analysis.

## Usage

``` r
prepare_CNV_SE(
  cnv_data,
  clinical_data,
  variability_threshold = 0.25,
  method = NULL,
  p_threshold = NULL,
  top_n_supervised = NULL
)
```

## Arguments

- cnv_data:

  A numeric matrix of CNV data where:

  - Rows represent genes/genomic regions

  - Columns represent samples

  - Values typically represent copy number states (e.g., -2, -1, 0,
    1, 2) or continuous copy number estimates

- clinical_data:

  A matrix or data frame of clinical data where:

  - Columns represent samples (matching `cnv_data` column names)

  - Rows represent clinical variables

  - Must include rows named `"overall_survival"` and `"status"` if using
    supervised feature selection

- variability_threshold:

  Numeric value between 0 and 1 specifying the proportion of most
  variable genes to retain when `method = "unsupervised"`. For example,
  0.25 retains the top 25\\ is not "unsupervised".

- method:

  Character string specifying the feature selection method. Options are:

  - `"supervised"`: Cox regression-based selection using survival
    outcomes

  - `"unsupervised"`: Variability-based selection (no outcome
    information)

  - `NULL`: No additional filtering beyond constant gene removal
    (default)

  Default is NULL.

- p_threshold:

  Numeric p-value threshold for supervised feature selection. Only genes
  with Cox regression p-values \<= `p_threshold` are retained. Required
  when `method = "supervised"`. Ignored otherwise. Default is NULL.

- top_n_supervised:

  Integer specifying the number of top genes to select if fewer than 100
  genes pass the p-value threshold in supervised selection. Acts as a
  fallback to ensure a minimum number of features. Default is NULL.

## Value

A SummarizedExperiment object containing:

- `assays`: List with one element named "counts" containing the filtered
  CNV matrix

- `metadata`: List with preprocessing information:

  - `omics`: "CNV"

  - `method`: The feature selection method used

## Details

The function performs the following processing steps:

1.  **Sample alignment**: Identifies and retains only samples present in
    both CNV and clinical data

2.  **Constant gene removal**: Removes genes with no variation across
    samples (all values identical)

3.  **Feature selection** (optional):

    - **Supervised**: Fits univariate Cox regression for each gene,
      retains genes with p-value \<= `p_threshold`. Uses parallel
      processing for speed. Falls back to top N genes if \< 100 pass
      threshold.

    - **Unsupervised**: Calculates variability as sum of absolute CNV
      values, retains top `variability_threshold` proportion of genes

4.  **SummarizedExperiment creation**: Packages filtered data with
    metadata

## See also

[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
for the return object class
[`coxph`](https://rdrr.io/pkg/survival/man/coxph.html) for Cox
regression in supervised selection
