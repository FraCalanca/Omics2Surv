# Build a C-index heatmap from a train/test MAE pair

Wrapper: harvests the metrics out of a `train` and `test` MAE object,
tidies/classifies them, and returns the faceted heatmap as a `ggplot`
object.

## Usage

``` r
c_heatmap(
  mae_train,
  mae_test,
  name = "Dataset",
  aftcoop_rho = NULL,
  strat_order = STRAT_ORDER,
  metric_order = METRIC_ORDER,
  std_methods = STD_METHODS,
  modality_map = MODALITY_MAP,
  branch_map = BRANCH_MAP,
  metric_map = METRIC_MAP,
  single_omics = SINGLE_OMICS,
  fill_limits = NULL,
  text_threshold = NULL,
  palette = "YlGnBu"
)
```

## Arguments

- mae_train, mae_test:

  MAE objects exposing `@metadata[["analysis"]]` for the training and
  test splits.

- name:

  Character scalar; the dataset label shown as the facet-column title.
  Defaults to `"Dataset"`.

- aftcoop_rho:

  Numeric scalar; the `beta_rho` value to select for AFTCoop results
  (e.g. `0.3`, `0.5`, `0.7`). When `NULL` (default) all available
  `beta_rho_*` sub-keys are kept as separate rows. When supplied, only
  the sub-key whose suffix equals `as.character(aftcoop_rho)` is
  retained; rows from other `rho` values are silently dropped.

- strat_order:

  Preferred top-to-bottom ordering of strategies. Defaults to Omics1,
  Omics2, Early, Late, Joint, Flexynesis. Missing entries are ignored;
  unlisted-but-present strategies are appended.

- metric_order:

  Preferred left-to-right ordering of metrics. Defaults to H-train,
  H-test, Uno. Handled the same way as `strat_order`.

- std_methods, modality_map, branch_map, metric_map, single_omics:

  Lookup tables used to classify and label paths. Override these if your
  path keys or desired labels differ.

- fill_limits:

  Numeric length-2 colour-scale limits. `NULL` (default) derives a
  padded range from the data.

- text_threshold:

  Value above which tile labels turn white for contrast. `NULL`
  (default) uses the midpoint of `fill_limits`.

- palette:

  ColorBrewer sequential palette name for the fill scale. Defaults to
  `"YlGnBu"`.

## Value

A `ggplot` object. Print it to display, or pass it to `ggsave()`
yourself if you want to export it.

## Details

The result is flexible: only the strategies, methods and metrics
actually present are drawn. `strat_order` and `metric_order` are treated
as *preferences* — entries not present are skipped, and anything present
but not listed is appended (with a message) rather than dropped. So a
run with a single omics and a single method (e.g. only `cox_lasso` on
Omics1) produces a one-row heatmap with no wasted space.
