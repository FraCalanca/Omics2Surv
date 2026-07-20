# Prepare Copy Number Variation Data as SummarizedExperiment

Processes CNV data by aligning samples with clinical data, removing
constant features, and optionally performing supervised, unsupervised,
or biological variable screening.

## Usage

``` r
prepare_CNV(
  cnv_data,
  clinical_data,
  method = NULL,
  supervised_method = c("pval", "absmg", "padj", "mgpval"),
  unsupervised_method = c("sparsity", "laplacian"),
  p_threshold = 0.05,
  max_zero_fraction = 0.75,
  feature_fraction = 0.25,
  top_n_features = NULL,
  gene_list = NULL,
  n_neighbors = NULL,
  on_constant = c("drop", "error", "keep")
)
```

## Arguments

- cnv_data:

  Numeric matrix of CNV data (genes/regions as rows, samples as
  columns).

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

  Character string specifying the unsupervised screening method:
  `"sparsity"` (fraction of zero entries) or `"laplacian"` (Laplacian
  score). Default is `"sparsity"`.

- p_threshold:

  Numeric p-value threshold for supervised selection. Used only when
  `supervised_method` is `"pval"` or `"padj"` and `top_n_features` is
  `NULL`. Default is 0.05.

- max_zero_fraction:

  Numeric proportion (0-1) for sparsity filtering: genes with a zero
  fraction above this value are removed. Used only when
  `unsupervised_method = "sparsity"` and `top_n_features` is `NULL`.
  Default is 0.75.

- feature_fraction:

  Numeric proportion (0-1) of features retained by Laplacian screening
  when `top_n_features` is `NULL`. Default is 0.25.

- top_n_features:

  Integer. If provided, overrides threshold/proportion filtering.
  **Mandatory** when `supervised_method` is `"absmg"` or `"mgpval"`.
  Default is `NULL`.

- gene_list:

  Character vector of feature names to retain. Used only when
  `method = "biological"`. Default is `NULL`.

- n_neighbors:

  Integer. Number of neighbors for Laplacian score computation
  (`unsupervised_method = "laplacian"`). Defaults to \\\approx
  \sqrt{n}\\ if `NULL`.

- on_constant:

  Character. After screening, any sample left constant (all-equal across
  the selected genes) is always repaired by adding back, for each such
  sample, the single original gene where it is most strongly altered, so
  that every retained sample has non-zero variance and downstream tools
  (e.g. flexynesis) do not drop it. Samples that are constant on the
  FULL input matrix carry no information and cannot be repaired;
  `on_constant` decides their fate: `"drop"` (default) removes them from
  the object, `"error"` stops and lists them, `"keep"` leaves them in
  (they stay constant).

## Value

A `SummarizedExperiment` object with:

- `assays$counts`: Filtered CNV matrix.

- `metadata$preprocessing`: List with `omics`, `method`, `sub_method`,
  and number of features retained.
