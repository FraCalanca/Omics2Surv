# Train Penalized Regression Survival Models on a MultiAssayExperiment Object

Fits penalized Cox or AFT regression models on one or more omics layers
from a training `MultiAssayExperiment`, supporting single-omic, early
integration, and late integration strategies. Clinical data (survival
outcome and optional covariates) are extracted directly from
`colData(mae_train)`. Results are stored in the returned MAE.

## Usage

``` r
train_pen_reg_models(
  mae_train,
  omic_list,
  integration = c("none", "early", "late"),
  model_method = c("cox_lasso", "cox_adaptive_lasso", "cox_elastic_net", "aft_pen"),
  late_method = c("RRA", "average"),
  clinical_vars = NULL,
  seed = 123,
  ...
)
```

## Arguments

- mae_train:

  A `MultiAssayExperiment` object containing the training data. Must
  have `colData` with `"overall_survival"` and `"status"` columns.

- omic_list:

  Character vector of assay names in `mae_train` to use (e.g. `"mRNA"`
  or `c("mRNA", "CNV")`).

- integration:

  Character string specifying the integration strategy:

  - `"none"`: fits one model per omic in `omic_list`.

  - `"early"`: concatenates all omics features into a single matrix and
    fits one model.

  - `"late"`: fits one model per omic independently, then combines
    predictions via `late_method`.

  Default is `"none"`.

- model_method:

  Character string specifying the penalized model: `"cox_lasso"`,
  `"cox_adaptive_lasso"`, `"cox_elastic_net"`, or `"aft_pen"`. Default
  is `"cox_lasso"`.

- late_method:

  Character string for late integration ensemble: `"RRA"` (Robust Rank
  Aggregation) or `"average"`. Used only when `integration = "late"`.
  Default is `"RRA"`.

- clinical_vars:

  Character vector of column names in `colData(mae_train)` to append as
  additional covariates to the design matrix. Default is `NULL`.

- seed:

  Integer seed for reproducibility. Default is 123.

- ...:

  Additional arguments passed to the underlying train function (e.g.
  `nfolds`, `optimal_lambda`, `alpha`).

## Value

`mae_train` with results stored in `metadata(mae_train)$analysis`:

- `"none"` / `"early"`: per-assay results at
  `[[model_method]][[assay_name]]`, containing `train` (model,
  predictions, selected features) and `metrics` (resubstitution Harrell
  C-index).

- `"late"`: per-omic results as above, plus ensemble results at
  `$late_integration[[model_method]][[tag]]` (where `tag` encodes the
  integrated omics), containing `train`, `metrics`, and `late_method`.

## See also

[`predict_pen_reg_models`](https://fracalanca.github.io/Omics2Surv/reference/predict_pen_reg_models.md),
[`build_MAE`](https://fracalanca.github.io/Omics2Surv/reference/build_MAE.md),
[`early_integration`](https://fracalanca.github.io/Omics2Surv/reference/early_integration.md),
[`late_integration`](https://fracalanca.github.io/Omics2Surv/reference/late_integration.md)
