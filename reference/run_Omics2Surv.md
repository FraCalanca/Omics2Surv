# Run the complete multi-omics survival analysis pipeline

Single-call wrapper that orchestrates the full workflow: optional data
retrieval from LinkedOmics, optional train/test split, per-omic feature
screening, MAE construction, model training (and prediction when a test
set exists), benchmarking, and downstream analyses (stratification,
Kaplan-Meier, feature heatmap, enrichment) on the best model.

## Usage

``` r
run_Omics2Surv(
  omics_list = NULL,
  clinical_data = NULL,
  cancer_type = NULL,
  linkedomics_omics = NULL,
  screening_specs,
  model_specs,
  train_fraction = 0.7,
  seed = 4820,
  best_by = c("uno_test", "harrell_test"),
  do_benchmark = TRUE,
  do_stratification = TRUE,
  strat_method = c("median", "survcutpoint", "user"),
  strat_threshold = NULL,
  do_km = TRUE,
  do_heatmap = TRUE,
  heatmap_annotations = NULL,
  do_enrichment = FALSE,
  enrichment_params = list(),
  aftcoop_rho = NULL,
  dataset_name = "Dataset",
  verbose = TRUE
)
```

## Arguments

- omics_list:

  Named list of raw omics matrices (features x samples), or `NULL` when
  using the LinkedOmics download mode. Names become assay names in the
  MAE (e.g. `list(mRNA = ..., CNV = ...)`).

- clinical_data:

  Matrix/data frame of clinical data (variables x samples), or `NULL` in
  LinkedOmics mode. Must contain rows `"overall_survival"` and
  `"status"`.

- cancer_type:

  Character. TCGA cancer-type string for LinkedOmics download (e.g.
  `"TCGA-LGG"`). Ignored when `omics_list` is supplied. Default `NULL`.

- linkedomics_omics:

  Named list mapping the assay names you want to the LinkedOmics
  data-type strings (e.g. `list(mRNA = "RNASeq (HiSeq, Gene Level)")`).
  Required in LinkedOmics mode. Default `NULL`.

- screening_specs:

  Named list (same names as `omics_list` / `linkedomics_omics`). Each
  element is a list describing how to screen that omic. See Details.

- model_specs:

  List of model requests. See Details.

- train_fraction:

  Numeric in (0, 1) for the train/test split, or `NULL` to skip
  splitting and train on the entire dataset. When `NULL`, no test set is
  built, no predictions are made, and `best_by` is forced to
  `"harrell_train"`. Default `0.7`.

- seed:

  Integer seed for the split and cross-validation. Default `4820`.

- best_by:

  Character. Metric used to pick the best model: `"uno_test"` or
  `"harrell_test"` when a test set exists; automatically forced to
  `"harrell_train"` when `train_fraction = NULL`. Default `"uno_test"`.

- do_benchmark:

  Logical. Build the benchmarking table? Default `TRUE`.

- do_stratification:

  Logical. Run patient stratification on the best model? Default `TRUE`.

- strat_method:

  Stratification cutpoint method: `"median"`, `"survcutpoint"`, or
  `"user"`. Default `"median"`.

- strat_threshold:

  Numeric threshold for `strat_method = "user"`. Default `NULL`.

- do_km:

  Logical. Draw Kaplan-Meier curves? Requires
  `do_stratification = TRUE`. Default `TRUE`.

- do_heatmap:

  Logical. Draw feature heatmap? Requires `do_stratification = TRUE`.
  Skipped automatically for late-integration ensembles. Default `TRUE`.

- heatmap_annotations:

  Character vector of `colData` column names to add as annotation tracks
  in the heatmap. Default `NULL`.

- do_enrichment:

  Logical. Run pathway enrichment on the best model's features? Skipped
  for late-integration ensembles. Default `FALSE`.

- enrichment_params:

  Named list of extra arguments forwarded to
  [`enrich_features`](https://fracalanca.github.io/Omics2Surv/reference/enrich_features.md).
  Default [`list()`](https://rdrr.io/r/base/list.html).

- aftcoop_rho:

  Numeric or `NULL`. For AFTCoop models: selects the `beta_rho_<rho>`
  variant in benchmark, heatmap, stratification, feature heatmap and
  enrichment. Auto-resolves when only one rho was trained; errors
  clearly when multiple are present. Default `NULL`.

- dataset_name:

  Character label used in
  [`c_heatmap`](https://fracalanca.github.io/Omics2Surv/reference/c_heatmap.md).
  Default `"Dataset"`.

- verbose:

  Logical. Print progress messages? Default `TRUE`.

## Value

A named list with:

- `mae_train`: MAE with all fitted models in its metadata.

- `mae_test`: MAE with predictions in its metadata, or `NULL` when
  `train_fraction = NULL`.

- `benchmark`: data frame from
  [`benchmark_models`](https://fracalanca.github.io/Omics2Surv/reference/benchmark_models.md),
  or `NULL`.

- `best_model`: named list with `model`, `modality`, `integration`,
  `metric`, and `value`.

- `stratification`: output of
  [`patient_stratification`](https://fracalanca.github.io/Omics2Surv/reference/patient_stratification.md),
  or `NULL`.

- `km_plot`: `ggsurvplot` object, or `NULL`.

- `heatmap`: pheatmap object, or `NULL`.

- `c_index_heatmap`: ggplot from
  [`c_heatmap`](https://fracalanca.github.io/Omics2Surv/reference/c_heatmap.md),
  or `NULL`.

- `enrichment`: output of
  [`enrich_features`](https://fracalanca.github.io/Omics2Surv/reference/enrich_features.md),
  or `NULL`.

## Details

### screening_specs

Each element must have a `type` field (`"RNA"`, `"methyl"`, or `"CNV"`)
and typically a `method` field (`"supervised"`, `"unsupervised"`,
`"biological"`, or `NULL`). All other fields are forwarded to
[`prepare_RNA()`](https://fracalanca.github.io/Omics2Surv/reference/prepare_RNA.md)
/
[`prepare_methyl()`](https://fracalanca.github.io/Omics2Surv/reference/prepare_methyl.md)
/
[`prepare_CNV()`](https://fracalanca.github.io/Omics2Surv/reference/prepare_CNV.md);
unrecognised fields are ignored (with a message), so a shared template
can be reused across omics.

### model_specs

Each element is a list with `type` (`"pen_reg"`, `"joint"`, or
`"flexynesis"`), `omic_list` (assay names to use), and type-specific
fields (`model_method`/`integration`/`late_method` for pen_reg;
`joint_method` for joint). Extra arguments go in a `params` sub-list
forwarded to the underlying train function.

### No-split mode (`train_fraction = NULL`)

The entire dataset is used for training. No test MAE is built, no
predictions are computed, and `best_by` is forced to `"harrell_train"`.
Downstream steps (stratification, KM, heatmap) run on `mae_train`
itself.

## Data input

Supply data in one of two ways:

- **User-supplied matrices** (default): pass `omics_list` and
  `clinical_data` directly. Any LinkedOmics-like matrices work here -
  from
  [`retrieve_linkedomics_data`](https://fracalanca.github.io/Omics2Surv/reference/retrieve_linkedomics_data.md),
  [`prepare_custom_data`](https://fracalanca.github.io/Omics2Surv/reference/prepare_custom_data.md),
  or any other source.

- **LinkedOmics download**: pass `cancer_type` and `linkedomics_omics`
  instead. The function calls
  [`retrieve_linkedomics_data`](https://fracalanca.github.io/Omics2Surv/reference/retrieve_linkedomics_data.md)
  automatically, so `omics_list` and `clinical_data` must be `NULL`.

## See also

[`retrieve_linkedomics_data`](https://fracalanca.github.io/Omics2Surv/reference/retrieve_linkedomics_data.md),
[`prepare_custom_data`](https://fracalanca.github.io/Omics2Surv/reference/prepare_custom_data.md),
[`split_train_test`](https://fracalanca.github.io/Omics2Surv/reference/split_train_test.md),
[`prepare_RNA`](https://fracalanca.github.io/Omics2Surv/reference/prepare_RNA.md),
[`build_MAE`](https://fracalanca.github.io/Omics2Surv/reference/build_MAE.md),
[`train_pen_reg_models`](https://fracalanca.github.io/Omics2Surv/reference/train_pen_reg_models.md),
[`benchmark_models`](https://fracalanca.github.io/Omics2Surv/reference/benchmark_models.md),
[`patient_stratification`](https://fracalanca.github.io/Omics2Surv/reference/patient_stratification.md),
[`plot_km`](https://fracalanca.github.io/Omics2Surv/reference/plot_km.md),
[`plot_heatmap`](https://fracalanca.github.io/Omics2Surv/reference/plot_heatmap.md),
[`enrich_features`](https://fracalanca.github.io/Omics2Surv/reference/enrich_features.md)
