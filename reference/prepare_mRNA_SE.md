# Prepare mRNA Expression Data as SummarizedExperiment

Processes mRNA expression data by aligning samples with clinical data,
removing genes with missing values or low expression, and optionally
performing supervised or unsupervised feature selection.

## Usage

``` r
prepare_mRNA_SE(
  RNA_data,
  clinical_data,
  q_threshold = 0.25,
  method = NULL,
  p_threshold = 0.1,
  hvg_threshold = 25
)
```

## Arguments

- RNA_data:

  Numeric matrix of mRNA expression data (genes as rows, samples as
  columns).

- clinical_data:

  Matrix or data frame of clinical data (variables as rows, samples as
  columns). Must include "overall_survival" and "status" rows for
  supervised selection.

- q_threshold:

  Numeric quantile threshold (0-1) for low expression filtering. Genes
  with expression at this quantile equal to 0 are removed. Default is
  0.25.

- method:

  Character string for feature selection: `"supervised"` (Cox
  regression-based), `"unsupervised"` (high variable genes), or `NULL`
  (no selection). Default is NULL.

- p_threshold:

  Numeric p-value threshold for supervised selection. Default is 0.1.

- hvg_threshold:

  Numeric percentile (0-100) for highly variable gene selection in
  unsupervised method. Default is 25 (top 25% most variable genes).

## Value

A SummarizedExperiment object with:

- `assays$counts`: Filtered mRNA expression matrix

- `metadata$preprocessing`: List with `omics = "mRNA"` and `method` used
