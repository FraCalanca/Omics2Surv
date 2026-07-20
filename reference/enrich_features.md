# Functional Enrichment of a Model's Features

Runs over-representation analysis (ORA) on the genes of a survival
model: Gene Ontology and KEGG via clusterProfiler, Reactome via
ReactomePA. Features are obtained exactly as in
[`plot_umap`](https://fracalanca.github.io/Omics2Surv/reference/plot_umap.md)
/
[`plot_heatmap`](https://fracalanca.github.io/Omics2Surv/reference/plot_heatmap.md)
—
[`selected_features()`](https://fracalanca.github.io/Omics2Surv/reference/selected_features.md)
for penalised models and
[`top_features()`](https://fracalanca.github.io/Omics2Surv/reference/top_features.md)
for importance-based models (`"blockForest"`, `"flexynesis"`) — or
supplied directly via `features`. Gene symbols are mapped to ENTREZ IDs
for KEGG and Reactome; GO is queried on the original key type.

## Usage

``` r
enrich_features(
  mae,
  model,
  modality = NULL,
  top_feat = NULL,
  features = NULL,
  rho = NULL,
  strip_prefix = NULL,
  databases = c("GO", "KEGG", "Reactome"),
  ont = "ALL",
  key_type = "SYMBOL",
  org_db = "org.Hs.eg.db",
  kegg_organism = "hsa",
  reactome_organism = "human",
  pvalue_cutoff = 0.05,
  qvalue_cutoff = 0.2,
  ...
)
```

## Arguments

- mae:

  A `MultiAssayExperiment` whose `metadata()$analysis` holds the fitted
  model (typically the training object).

- model:

  Character string naming the model (e.g. `"cox_elastic_net"`,
  `"blockForest"`, `"flexynesis"`).

- modality:

  Character string naming the omics modality used to locate a penalised
  model's selected features (ignored for `"blockForest"` /
  `"flexynesis"`). Default is `NULL`: the single per-run integration tag
  (`"early_<omics>"`) is resolved automatically; pass the tag (or an
  omic name) explicitly to disambiguate.

- top_feat:

  Integer number of top features. Used *only* for `"blockForest"` /
  `"flexynesis"`. If `NULL` (default),
  [`top_features`](https://fracalanca.github.io/Omics2Surv/reference/top_features.md)'s
  own default is used.

- features:

  Optional character vector of gene identifiers, overriding the
  automatic retrieval from `model`.

- rho:

  Numeric. For `model = "AFTCoop"` only: selects the `beta_rho_<rho>`
  feature set. If `NULL` (default) and only one rho was trained, it is
  resolved automatically.

- strip_prefix:

  Optional regular expression removed from the start of each feature
  name before mapping (e.g. `"^(mRNA|meth)_"` when assay rows are
  layer-prefixed). Default is `NULL`.

- databases:

  Which analyses to run; any of `"GO"`, `"KEGG"`, `"Reactome"`. Default
  is all three.

- ont:

  GO ontology: `"BP"`, `"MF"`, `"CC"` or `"ALL"`. Default is `"ALL"`.

- key_type:

  Key type of the input features for GO and ID mapping. Default is
  `"SYMBOL"`.

- org_db:

  Name of the organism annotation package. Default is `"org.Hs.eg.db"`.

- kegg_organism:

  KEGG organism code. Default is `"hsa"`.

- reactome_organism:

  ReactomePA organism. Default is `"human"`.

- pvalue_cutoff, qvalue_cutoff:

  Significance cut-offs passed to each enrichment call. Defaults are
  0.05 and 0.2.

- ...:

  Additional arguments forwarded to the enrichment functions.

## Value

A named `list` of `enrichResult` objects (one per requested database),
plus `genes` (the cleaned input), `entrez` (the mapped IDs) and
`universe` (the background used, or `NULL`). KEGG and Reactome results
are made `readable` (ENTREZ back to symbols).

## See also

[`plot_enrichment`](https://fracalanca.github.io/Omics2Surv/reference/plot_enrichment.md),
[`selected_features`](https://fracalanca.github.io/Omics2Surv/reference/selected_features.md),
[`top_features`](https://fracalanca.github.io/Omics2Surv/reference/top_features.md)
