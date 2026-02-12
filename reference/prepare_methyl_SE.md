# Prepare Methylation Data as SummarizedExperiment

Processes methylation data by aligning samples with clinical data,
removing incomplete cases, and optionally performing supervised or
unsupervised feature selection.

## Usage

``` r
prepare_methyl_SE(
  methyl_data,
  clinical_data,
  system_type = "unix",
  method = NULL,
  p_threshold = 0.1,
  iqr_threshold = 0.25
)
```

## Arguments

- methyl_data:

  Numeric matrix of methylation data (genes/CpG sites as rows, samples
  as columns).

- clinical_data:

  Matrix or data frame of clinical data (variables as rows, samples as
  columns). Must include "overall_survival" and "status" rows for
  supervised selection.

- system_type:

  Character string specifying OS type ("unix" or other). Currently
  unused; parallel backend auto-detected from platform.

- method:

  Character string for feature selection: `"supervised"` (Cox
  regression-based), `"unsupervised"` (IQR-based), or `NULL` (no
  selection). Default is NULL.

- p_threshold:

  Numeric p-value threshold for supervised selection. Default is 0.1.

- iqr_threshold:

  Numeric proportion (0-1) of top IQR genes to retain for unsupervised
  selection. Default is 0.25.

## Value

A SummarizedExperiment object with:

- `assays$counts`: Filtered methylation matrix

- `metadata$preprocessing`: List with `omics = "methylation"` and
  `method` used
