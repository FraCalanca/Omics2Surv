# Prepare RNA Expression Data as SummarizedExperiment

Processes RNA expression data (e.g., mRNA or miRNA) by aligning samples
with clinical data, removing features with low expression, and
optionally performing supervised, unsupervised, or biological variable
screening.

## Usage

``` r
prepare_RNA(
  RNA_data,
  clinical_data,
  q_threshold = 0.25,
  method = NULL,
  supervised_method = c("pval", "absmg", "padj", "mgpval"),
  unsupervised_method = c("A", "MAD", "laplacian"),
  p_threshold = 0.05,
  hvg_threshold = 25,
  top_n_features = NULL,
  gene_list = NULL,
  n_neighbors = NULL
)
```

## Arguments

- RNA_data:

  Numeric matrix of RNA expression data (features as rows, samples as
  columns). Features may represent mRNAs, miRNAs, or other
  transcriptomic measurements.

- clinical_data:

  Matrix or data frame of clinical data (variables as rows, samples as
  columns). Must include "overall_survival" and "status" rows for
  supervised selection.

- q_threshold:

  Numeric quantile threshold (0-1) for low-expression filtering.
  Features with expression at this quantile equal to 0 are removed.
  Default is 0.25.

- method:

  Character string for variable screening: `"supervised"`,
  `"unsupervised"`, `"biological"`, or `NULL` (no selection). Default is
  `NULL`.

- supervised_method:

  Character string specifying the supervised screening criterion, passed
  to
  [`supervised_screening`](https://fracalanca.github.io/Omics2Surv/reference/supervised_screening.md):

  - `"pval"`: p-value of the marginal Cox test. Supports `p_threshold`
    and `top_n_features`.

  - `"padj"`: Benjamini-Hochberg adjusted p-value. Supports
    `p_threshold` and `top_n_features`.

  - `"absmg"`: absolute marginal Cox coefficient \\\|\hat\beta_j\|\\.
    Requires `top_n_features`.

  - `"mgpval"`: combined score \\-\log\_{10}(p)\times\|\hat\beta_j\|\\.
    Requires `top_n_features`.

  Default is `"pval"`.

- unsupervised_method:

  Character string specifying the unsupervised screening method:

  - `"A"`: coefficient of variation.

  - `"MAD"`: median absolute deviation.

  - `"laplacian"`: Laplacian score.

  Default is `"A"`.

- p_threshold:

  Numeric p-value threshold for supervised selection. Used only when
  `supervised_method` is `"pval"` or `"padj"` and `top_n_features` is
  `NULL`. Default is 0.05.

- hvg_threshold:

  Numeric percentile (0-100) for highly variable feature selection in
  unsupervised methods. Used only when `top_n_features` is `NULL`.
  Default is 25 (top 25% most variable features).

- top_n_features:

  Integer. If provided, overrides threshold-based filtering. Mandatory
  when `supervised_method` is `"absmg"` or `"mgpval"`. Default is
  `NULL`.

- gene_list:

  Character vector of feature names to retain. Used only when
  `method = "biological"`. Features may correspond to genes, miRNAs, or
  other transcriptomic entities. Default is `NULL`.

- n_neighbors:

  Integer. Number of neighbors for Laplacian score computation
  (`unsupervised_method = "laplacian"`). Recommended to tune based on
  sample size, approximately \\\sqrt{n}\\. Default is `NULL`.

## Value

A `SummarizedExperiment` object containing:

- `assays$counts`: Filtered RNA expression matrix.

- `metadata$preprocessing`: A list containing:

  - `omics`: Omics type. Always `"RNA"` (covers mRNA, miRNA, and other
    transcriptomic measurements).

  - `method`: Screening method applied.

  - `sub_method`: Specific screening criterion used.

  - `n_features`: Number of retained features.
