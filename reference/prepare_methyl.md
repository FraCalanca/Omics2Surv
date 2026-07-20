# Prepare Methylation Data as SummarizedExperiment

Processes methylation data by aligning samples with clinical data, and
optionally performing supervised, unsupervised, or biological feature
screening.

## Usage

``` r
prepare_methyl(
  methyl_data,
  clinical_data,
  method = NULL,
  supervised_method = c("pval", "absmg", "padj", "mgpval"),
  unsupervised_method = c("IQR", "MAD", "laplacian"),
  p_threshold = 0.05,
  hvg_threshold = 25,
  top_n_features = NULL,
  gene_list = NULL,
  n_neighbors = NULL
)
```

## Arguments

- methyl_data:

  Numeric matrix or data.frame of methylation data (features as rows,
  samples as columns).

- clinical_data:

  Matrix or data frame of clinical data (variables as rows, samples as
  columns). Must include `"overall_survival"` and `"status"` rows for
  supervised selection.

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
    **Requires `top_n_features`**.

  - `"mgpval"`: combined score \\-\log\_{10}(p) \times
    \|\hat\beta_j\|\\. **Requires `top_n_features`**.

  Default is `"pval"`.

- unsupervised_method:

  Character string specifying the unsupervised screening method: `"IQR"`
  (interquartile range), `"MAD"` (median absolute deviation), or
  `"laplacian"` (Laplacian score). Default is `"IQR"`.

- p_threshold:

  Numeric p-value threshold for supervised selection. Used only when
  `supervised_method` is `"pval"` or `"padj"` and `top_n_features` is
  `NULL`. Default is 0.05.

- hvg_threshold:

  Numeric percentile (0-100) of top features to retain in unsupervised
  methods. Used only when `top_n_features` is `NULL`. Default is 25 (top
  25%).

- top_n_features:

  Integer. If provided, overrides threshold-based filtering.
  **Mandatory** when `supervised_method` is `"absmg"` or `"mgpval"`.
  Default is `NULL`.

- gene_list:

  Character vector of gene names to retain. Used only when
  `method = "biological"`. Default is `NULL`.

- n_neighbors:

  Integer. Number of neighbors for Laplacian score computation
  (`unsupervised_method = "laplacian"`). Defaults to \\\approx
  \sqrt{n}\\ if `NULL`.

## Value

A `SummarizedExperiment` object with:

- `assays$counts`: Filtered methylation matrix.

- `metadata$preprocessing`: List with `omics`, `method`, `sub_method`,
  and number of features retained.
