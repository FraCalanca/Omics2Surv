# =============================================================================
# Internal helpers
# =============================================================================

# Null-coalescing operator.
`%||%` <- function(x, y) if (is.null(x)) y else x

# Joint model names as stored by train_j_models / benchmark_models.
#' @keywords internal
#' @noRd
.is_joint_model <- function(model) {
  model %in% c("blockForest", "CoopCox", "AFTCoop") ||
    startsWith(model, "AFTCoop_")
}

# Container in metadata(mae)$analysis that holds a result, given the
# integration label reported by benchmark_models().
#' @keywords internal
#' @noRd
.result_container <- function(integration) {
  switch(integration,
         late  = "late_integration",
         joint = "joint",
         NULL)
}

# benchmark_models() labels AFTCoop rows as "AFTCoop_beta_rho_<v>" when
# aftcoop_rho is NULL; map back to "AFTCoop" so downstream readers find it.
#' @keywords internal
#' @noRd
.normalise_model_name <- function(model) {
  if (startsWith(model, "AFTCoop_")) "AFTCoop" else model
}


#' Run the complete multi-omics survival analysis pipeline
#'
#' Single-call wrapper that orchestrates the full workflow: optional data
#' retrieval from LinkedOmics, optional train/test split, per-omic feature
#' screening, MAE construction, model training (and prediction when a test set
#' exists), benchmarking, and downstream analyses (stratification,
#' Kaplan-Meier, feature heatmap, enrichment) on the best model.
#'
#' @section Data input:
#' Supply data in one of two ways:
#' \itemize{
#'   \item \strong{User-supplied matrices} (default): pass \code{omics_list}
#'     and \code{clinical_data} directly. Any LinkedOmics-like matrices work
#'     here - from \code{\link{retrieve_linkedomics_data}},
#'     \code{\link{prepare_custom_data}}, or any other source.
#'   \item \strong{LinkedOmics download}: pass \code{cancer_type} and
#'     \code{linkedomics_omics} instead. The function calls
#'     \code{\link{retrieve_linkedomics_data}} automatically, so
#'     \code{omics_list} and \code{clinical_data} must be \code{NULL}.
#' }
#'
#' @param omics_list Named list of raw omics matrices (features x samples), or
#'   \code{NULL} when using the LinkedOmics download mode. Names become assay
#'   names in the MAE (e.g. \code{list(mRNA = ..., CNV = ...)}).
#' @param clinical_data Matrix/data frame of clinical data (variables x
#'   samples), or \code{NULL} in LinkedOmics mode. Must contain rows
#'   \code{"overall_survival"} and \code{"status"}.
#' @param cancer_type Character. TCGA cancer-type string for LinkedOmics
#'   download (e.g. \code{"TCGA-LGG"}). Ignored when \code{omics_list} is
#'   supplied. Default \code{NULL}.
#' @param linkedomics_omics Named list mapping the assay names you want to the
#'   LinkedOmics data-type strings (e.g.
#'   \code{list(mRNA = "RNASeq (HiSeq, Gene Level)")}). Required in LinkedOmics
#'   mode. Default \code{NULL}.
#' @param screening_specs Named list (same names as \code{omics_list} /
#'   \code{linkedomics_omics}). Each element is a list describing how to screen
#'   that omic. See Details.
#' @param model_specs List of model requests. See Details.
#' @param train_fraction Numeric in (0, 1) for the train/test split, or
#'   \code{NULL} to skip splitting and train on the entire dataset. When
#'   \code{NULL}, no test set is built, no predictions are made, and
#'   \code{best_by} is forced to \code{"harrell_train"}. Default \code{0.7}.
#' @param seed Integer seed for the split and cross-validation. Default
#'   \code{4820}.
#' @param best_by Character. Metric used to pick the best model:
#'   \code{"uno_test"} or \code{"harrell_test"} when a test set exists;
#'   automatically forced to \code{"harrell_train"} when
#'   \code{train_fraction = NULL}. Default \code{"uno_test"}.
#' @param do_benchmark Logical. Build the benchmarking table? Default
#'   \code{TRUE}.
#' @param do_stratification Logical. Run patient stratification on the best
#'   model? Default \code{TRUE}.
#' @param strat_method Stratification cutpoint method:
#'   \code{"median"}, \code{"survcutpoint"}, or \code{"user"}.
#'   Default \code{"median"}.
#' @param strat_threshold Numeric threshold for \code{strat_method = "user"}.
#'   Default \code{NULL}.
#' @param do_km Logical. Draw Kaplan-Meier curves? Requires
#'   \code{do_stratification = TRUE}. Default \code{TRUE}.
#' @param do_heatmap Logical. Draw feature heatmap? Requires
#'   \code{do_stratification = TRUE}. Skipped automatically for
#'   late-integration ensembles. Default \code{TRUE}.
#' @param heatmap_annotations Character vector of \code{colData} column names
#'   to add as annotation tracks in the heatmap. Default \code{NULL}.
#' @param do_enrichment Logical. Run pathway enrichment on the best model's
#'   features? Skipped for late-integration ensembles. Default \code{FALSE}.
#' @param enrichment_params Named list of extra arguments forwarded to
#'   \code{\link{enrich_features}}. Default \code{list()}.
#' @param aftcoop_rho Numeric or \code{NULL}. For AFTCoop models: selects the
#'   \code{beta_rho_<rho>} variant in benchmark, heatmap, stratification,
#'   feature heatmap and enrichment. Auto-resolves when only one rho was
#'   trained; errors clearly when multiple are present. Default \code{NULL}.
#' @param dataset_name Character label used in \code{\link{c_heatmap}}.
#'   Default \code{"Dataset"}.
#' @param verbose Logical. Print progress messages? Default \code{TRUE}.
#'
#' @return A named list with:
#'   \itemize{
#'     \item \code{mae_train}: MAE with all fitted models in its metadata.
#'     \item \code{mae_test}: MAE with predictions in its metadata, or
#'       \code{NULL} when \code{train_fraction = NULL}.
#'     \item \code{benchmark}: data frame from \code{\link{benchmark_models}},
#'       or \code{NULL}.
#'     \item \code{best_model}: named list with \code{model}, \code{modality},
#'       \code{integration}, \code{metric}, and \code{value}.
#'     \item \code{stratification}: output of
#'       \code{\link{patient_stratification}}, or \code{NULL}.
#'     \item \code{km_plot}: \code{ggsurvplot} object, or \code{NULL}.
#'     \item \code{heatmap}: pheatmap object, or \code{NULL}.
#'     \item \code{c_index_heatmap}: ggplot from \code{\link{c_heatmap}},
#'       or \code{NULL}.
#'     \item \code{enrichment}: output of \code{\link{enrich_features}},
#'       or \code{NULL}.
#'   }
#'
#' @details
#' ## screening_specs
#' Each element must have a \code{type} field (\code{"RNA"}, \code{"methyl"},
#' or \code{"CNV"}) and typically a \code{method} field (\code{"supervised"},
#' \code{"unsupervised"}, \code{"biological"}, or \code{NULL}). All other
#' fields are forwarded to \code{prepare_RNA()} / \code{prepare_methyl()} /
#' \code{prepare_CNV()}; unrecognised fields are ignored (with a message), so
#' a shared template can be reused across omics.
#'
#' ## model_specs
#' Each element is a list with \code{type} (\code{"pen_reg"}, \code{"joint"},
#' or \code{"flexynesis"}), \code{omic_list} (assay names to use), and
#' type-specific fields (\code{model_method}/\code{integration}/\code{late_method}
#' for pen_reg; \code{joint_method} for joint). Extra arguments go in a
#' \code{params} sub-list forwarded to the underlying train function.
#'
#' ## No-split mode (\code{train_fraction = NULL})
#' The entire dataset is used for training. No test MAE is built, no
#' predictions are computed, and \code{best_by} is forced to
#' \code{"harrell_train"}. Downstream steps (stratification, KM, heatmap)
#' run on \code{mae_train} itself.
#'
#' @seealso \code{\link{retrieve_linkedomics_data}},
#'   \code{\link{prepare_custom_data}}, \code{\link{split_train_test}},
#'   \code{\link{prepare_RNA}}, \code{\link{build_MAE}},
#'   \code{\link{train_pen_reg_models}}, \code{\link{benchmark_models}},
#'   \code{\link{patient_stratification}}, \code{\link{plot_km}},
#'   \code{\link{plot_heatmap}}, \code{\link{enrich_features}}
#'
#' @export
run_Omics2Surv <- function(
    omics_list           = NULL,
    clinical_data        = NULL,
    cancer_type          = NULL,
    linkedomics_omics    = NULL,
    screening_specs,
    model_specs,
    train_fraction       = 0.7,
    seed                 = 4820,
    best_by              = c("uno_test", "harrell_test"),
    do_benchmark         = TRUE,
    do_stratification    = TRUE,
    strat_method         = c("median", "survcutpoint", "user"),
    strat_threshold      = NULL,
    do_km                = TRUE,
    do_heatmap           = TRUE,
    heatmap_annotations  = NULL,
    do_enrichment        = FALSE,
    enrichment_params    = list(),
    aftcoop_rho          = NULL,
    dataset_name         = "Dataset",
    verbose              = TRUE
) {

  strat_method <- match.arg(strat_method)
  best_by      <- match.arg(best_by)
  .msg <- function(...) if (verbose) message(...)

  # No-split mode: override best_by
  do_split <- !is.null(train_fraction)
  if (!do_split) {
    if (best_by %in% c("uno_test", "harrell_test")) {
      .msg("  (train_fraction = NULL: no test set - best_by forced to 'harrell_train')")
      best_by <- "harrell_train"
    }
  }


  # ---------------------------------------------------------------------------
  # 1. Data acquisition
  # ---------------------------------------------------------------------------
  .msg("== [1/7] Acquiring data ==")

  # Determine mode: user-supplied vs LinkedOmics download
  use_linkedomics <- !is.null(cancer_type) && !is.null(linkedomics_omics)
  use_supplied    <- !is.null(omics_list)  && !is.null(clinical_data)

  if (use_linkedomics && use_supplied)
    stop("Supply EITHER (omics_list + clinical_data) OR (cancer_type + ",
         "linkedomics_omics), not both.")
  if (!use_linkedomics && !use_supplied)
    stop("You must supply data: either (omics_list + clinical_data) or ",
         "(cancer_type + linkedomics_omics).")

  if (use_linkedomics) {
    if (!is.list(linkedomics_omics) || is.null(names(linkedomics_omics)))
      stop("'linkedomics_omics' must be a named list: ",
           "list(assay_name = 'LinkedOmics data type string', ...).")

    .msg("  Downloading from LinkedOmics: ", cancer_type)
    omics_list <- lapply(linkedomics_omics, function(data_type) {
      .msg("    Fetching: ", data_type)
      retrieve_linkedomics_data(cancer_type, data_type)
    })
    names(omics_list) <- names(linkedomics_omics)

    .msg("    Fetching: Clinical")
    clinical_data <- retrieve_linkedomics_data(cancer_type, "Clinical")

  } else {
    .msg("  Using user-supplied data (",
         length(omics_list), " omic(s), ",
         ncol(clinical_data), " samples in clinical).")
  }


  # ---------------------------------------------------------------------------
  # 2. Validate inputs
  # ---------------------------------------------------------------------------
  .msg("== [2/7] Validating inputs ==")

  if (!is.list(omics_list) || is.null(names(omics_list)))
    stop("'omics_list' must be a named list.")
  if (!is.list(screening_specs) || is.null(names(screening_specs)))
    stop("'screening_specs' must be a named list.")

  missing_specs <- setdiff(names(omics_list), names(screening_specs))
  if (length(missing_specs))
    stop("'screening_specs' is missing entries for: ",
         paste(missing_specs, collapse = ", "))

  for (nm in names(omics_list)) {
    sp <- screening_specs[[nm]]
    if (is.null(sp$type) || !sp$type %in% c("RNA", "methyl", "CNV"))
      stop("screening_specs[['", nm, "']]$type must be 'RNA', 'methyl', or 'CNV'.")
  }

  if (!is.list(model_specs) || length(model_specs) == 0L)
    stop("'model_specs' must be a non-empty list.")

  for (i in seq_along(model_specs)) {
    spec <- model_specs[[i]]
    if (is.null(spec$type) || !spec$type %in% c("pen_reg", "joint", "flexynesis"))
      stop("model_specs[[", i, "]]$type must be 'pen_reg', 'joint', or 'flexynesis'.")
    if (is.null(spec$omic_list))
      stop("model_specs[[", i, "]] is missing $omic_list.")
    bad_omics <- setdiff(spec$omic_list, names(omics_list))
    if (length(bad_omics))
      stop("model_specs[[", i, "]]$omic_list references unknown omics: ",
           paste(bad_omics, collapse = ", "))
  }


  # ---------------------------------------------------------------------------
  # 3. Train / test split  (optional)
  # ---------------------------------------------------------------------------
  prep_fun_map <- list(RNA = prepare_RNA, methyl = prepare_methyl, CNV = prepare_CNV)
  data_arg_map <- c(RNA = "RNA_data", methyl = "methyl_data", CNV = "cnv_data")

  mae_test       <- NULL
  train_se_list  <- vector("list", length(omics_list))
  names(train_se_list) <- names(omics_list)

  if (do_split) {

    .msg("== [3/7] Splitting train / test (fraction = ", train_fraction, ") ==")

    split <- split_train_test(
      omics_list     = omics_list,
      clinical_data  = clinical_data,
      train_fraction = train_fraction,
      seed           = seed
    )

    # -- Screening: train SEs + aligned test SEs ------------------------------
    .msg("== Feature screening ==")
    test_se_list        <- train_se_list

    for (nm in names(omics_list)) {
      spec      <- screening_specs[[nm]]
      omic_type <- spec$type
      .msg("  Screening '", nm, "' (type = ", omic_type,
           ", method = ", spec$method %||% "NULL", ")")

      prep_fun <- prep_fun_map[[omic_type]]
      data_arg <- data_arg_map[[omic_type]]

      prep_args               <- spec[setdiff(names(spec), "type")]
      prep_args$clinical_data <- split$train_clinical
      prep_args[[data_arg]]   <- split$train_list[[nm]]

      unused <- setdiff(names(prep_args), names(formals(prep_fun)))
      if (length(unused))
        .msg("    (ignoring: ", paste(unused, collapse = ", "), ")")
      prep_args <- prep_args[names(prep_args) %in% names(formals(prep_fun))]

      train_se            <- do.call(prep_fun, prep_args)
      train_se_list[[nm]] <- train_se

      test_se_list[[nm]]  <- prepare_test_SE(
        test_data = split$test_list[[nm]],
        train_se  = train_se
      )
    }

    # -- Build both MAEs -------------------------------------------------------
    .msg("== Building MAE objects ==")
    mae_train <- build_MAE(omics_list = train_se_list,
                           clinical_data = split$train_clinical)
    mae_test  <- build_MAE(omics_list = test_se_list,
                           clinical_data = split$test_clinical)

  } else {

    .msg("== [3/7] No split - using full dataset for training ==")

    # -- Screening: full dataset ----------------------------------------------
    .msg("== Feature screening ==")

    for (nm in names(omics_list)) {
      spec      <- screening_specs[[nm]]
      omic_type <- spec$type
      .msg("  Screening '", nm, "' (type = ", omic_type,
           ", method = ", spec$method %||% "NULL", ")")

      prep_fun <- prep_fun_map[[omic_type]]
      data_arg <- data_arg_map[[omic_type]]

      prep_args               <- spec[setdiff(names(spec), "type")]
      prep_args$clinical_data <- clinical_data
      prep_args[[data_arg]]   <- omics_list[[nm]]

      unused <- setdiff(names(prep_args), names(formals(prep_fun)))
      if (length(unused))
        .msg("    (ignoring: ", paste(unused, collapse = ", "), ")")
      prep_args <- prep_args[names(prep_args) %in% names(formals(prep_fun))]

      train_se_list[[nm]] <- do.call(prep_fun, prep_args)
    }

    # -- Build training MAE only -----------------------------------------------
    .msg("== Building MAE object ==")
    mae_train <- build_MAE(omics_list = train_se_list,
                           clinical_data = clinical_data)
  }


  # ---------------------------------------------------------------------------
  # 4-5. Train models  (+  predict on test when split exists)
  # ---------------------------------------------------------------------------
  .msg("== [4/7] Training", if (do_split) " and predicting" else "", " models ==")

  for (i in seq_along(model_specs)) {
    spec   <- model_specs[[i]]
    params <- spec$params %||% list()

    if (spec$type == "pen_reg") {
      .msg("  [pen_reg] method = ", spec$model_method %||% "cox_lasso",
           ", integration = ", spec$integration %||% "none",
           ", omics = ", paste(spec$omic_list, collapse = "+"))

      train_args <- c(
        list(mae_train     = mae_train,
             omic_list     = spec$omic_list,
             integration   = spec$integration  %||% "none",
             model_method  = spec$model_method %||% "cox_lasso",
             late_method   = spec$late_method  %||% "RRA",
             clinical_vars = spec$clinical_vars,
             seed          = seed),
        params
      )
      mae_train <- do.call(train_pen_reg_models, train_args)

      if (do_split) {
        mae_test <- predict_pen_reg_models(
          mae_pred      = mae_test,
          mae_train     = mae_train,
          omic_list     = spec$omic_list,
          integration   = spec$integration  %||% "none",
          model_method  = spec$model_method %||% "cox_lasso",
          clinical_vars = spec$clinical_vars
        )
      }

    } else if (spec$type == "joint") {
      if (!exists("train_j_models", mode = "function") ||
          !exists("predict_j_models", mode = "function"))
        stop("model_specs[[", i, "]] requests a joint model but ",
             "train_j_models()/predict_j_models() are not available.")

      .msg("  [joint] method = ", spec$joint_method %||% "blockForest",
           ", omics = ", paste(spec$omic_list, collapse = "+"))

      train_args <- c(
        list(mae_train    = mae_train,
             omic_list    = spec$omic_list,
             joint_method = spec$joint_method %||% "blockForest",
             seed         = seed),
        params
      )
      mae_train <- do.call(train_j_models, train_args)

      if (do_split) {
        pred_args <- c(
          list(mae_pred     = mae_test,
               mae_train    = mae_train,
               omic_list    = spec$omic_list,
               joint_method = spec$joint_method %||% "blockForest"),
          params[intersect(names(params), "threads")]
        )
        mae_test <- do.call(predict_j_models, pred_args)
      }

    } else if (spec$type == "flexynesis") {
      if (!exists("train_flexynesis_models", mode = "function") ||
          !exists("predict_flexynesis_models", mode = "function"))
        stop("model_specs[[", i, "]] requests flexynesis but ",
             "train_flexynesis_models()/predict_flexynesis_models() are not available.")

      .msg("  [flexynesis] omics = ", paste(spec$omic_list, collapse = "+"))

      mae_train <- do.call(train_flexynesis_models,
                           c(list(mae_train  = mae_train,
                                  omic_list  = spec$omic_list),
                             params))
      if (do_split) {
        mae_test <- do.call(predict_flexynesis_models,
                            c(list(mae_pred   = mae_test,
                                   mae_train  = mae_train,
                                   omic_list  = spec$omic_list),
                              params))
      }
    }
  }


  # ---------------------------------------------------------------------------
  # 5. Benchmarking + best-model selection
  # ---------------------------------------------------------------------------
  benchmark_tbl   <- NULL
  c_index_heatmap <- NULL
  best_model_info <- NULL

  if (do_benchmark) {
    .msg("== [5/7] Benchmarking (best_by = '", best_by, "') ==")

    benchmark_tbl <- tryCatch(
      if (do_split)
        benchmark_models(mae_train, mae_test, aftcoop_rho = aftcoop_rho)
      else
        benchmark_models(mae_train, mae_train, aftcoop_rho = aftcoop_rho),
      error = function(e) {
        warning("benchmark_models() failed: ", conditionMessage(e)); NULL
      }
    )

    # C-index heatmap only when a real test set exists
    if (do_split) {
      c_index_heatmap <- tryCatch(
        c_heatmap(mae_train, mae_test, name = dataset_name,
                    aftcoop_rho = aftcoop_rho),
        error = function(e) {
          warning("c_heatmap() failed: ", conditionMessage(e)); NULL
        }
      )
    }

    if (!is.null(benchmark_tbl) && nrow(benchmark_tbl) > 0L) {
      metric_vals <- benchmark_tbl[[best_by]]
      valid       <- which(!is.na(metric_vals))
      if (length(valid)) {
        best_row <- benchmark_tbl[valid[which.max(metric_vals[valid])], ,
                                  drop = FALSE]
        best_model_info <- list(
          model       = best_row$model,
          modality    = best_row$dataset,
          integration = best_row$integration,
          metric      = best_by,
          value       = best_row[[best_by]]
        )
        .msg("  Best model: ", best_model_info$model,
             " / ", best_model_info$modality,
             "  (", best_by, " = ", round(best_model_info$value, 3), ")")
      } else {
        warning("No non-missing '", best_by, "' values; downstream steps skipped.")
      }
    } else {
      warning("Empty benchmark table; downstream steps skipped.")
    }
  }


  # ---------------------------------------------------------------------------
  # 6. Downstream: stratification, KM, feature heatmap, enrichment
  # ---------------------------------------------------------------------------
  stratification <- NULL
  km_plot        <- NULL
  heatmap_plot   <- NULL
  enrichment_res <- NULL

  # The MAE used for downstream viz: test when available, else train.
  mae_viz <- if (do_split) mae_test else mae_train

  if (!is.null(best_model_info) &&
      (do_stratification || do_km || do_heatmap || do_enrichment)) {

    model       <- .normalise_model_name(best_model_info$model)
    modality    <- best_model_info$modality
    integration <- best_model_info$integration
    container   <- .result_container(integration)
    feature_ok  <- !identical(integration, "late")

    .msg("== [6/7] Downstream analyses (best: ", model, " / ", modality, ") ==")

    # -- 6a. Patient stratification -------------------------------------------
    if (do_stratification) {
      .msg("  Patient stratification ...")
      strat_call <- list(
        mae_train  = mae_train,
        mae_test   = if (do_split) mae_test else NULL,
        model      = model,
        modality   = modality,
        container  = container,
        rho        = aftcoop_rho,
        method     = strat_method,
        threshold  = strat_threshold
      )
      stratification <- tryCatch(
        do.call(patient_stratification, strat_call),
        error = function(e) {
          warning("patient_stratification() failed: ", conditionMessage(e))
          NULL
        }
      )
    }

    # -- 6b. Kaplan-Meier plot ------------------------------------------------
    if (do_km && !is.null(stratification)) {
      .msg("  Kaplan-Meier plot ...")
      km_plot <- tryCatch(
        plot_km(stratification),
        error = function(e) { warning("plot_km() failed: ", conditionMessage(e)); NULL }
      )
    }

    # -- 6c. Feature heatmap --------------------------------------------------
    if (do_heatmap && !is.null(stratification) && feature_ok) {
      .msg("  Feature heatmap ...")
      heatmap_plot <- tryCatch(
        plot_heatmap(
          mae         = mae_viz,
          x           = stratification,
          model       = model,
          modality    = modality,
          rho         = aftcoop_rho,
          mae_train   = mae_train,
          annotations = heatmap_annotations
        ),
        error = function(e) {
          warning("plot_heatmap() failed: ", conditionMessage(e)); NULL
        }
      )
    } else if (do_heatmap && !feature_ok) {
      .msg("  Feature heatmap skipped (late-integration ensemble).")
    }

    # -- 6d. Enrichment analysis ----------------------------------------------
    if (do_enrichment && feature_ok) {
      .msg("  Enrichment analysis ...")
      enrich_call <- c(
        list(mae = mae_train, model = model, modality = modality,
             rho = aftcoop_rho),
        enrichment_params
      )
      enrichment_res <- tryCatch(
        do.call(enrich_features, enrich_call),
        error = function(e) {
          warning("enrich_features() failed: ", conditionMessage(e)); NULL
        }
      )
    } else if (do_enrichment && !feature_ok) {
      .msg("  Enrichment skipped (late-integration ensemble).")
    }
  }

  .msg("== [7/7] Pipeline complete ==")

  list(
    mae_train       = mae_train,
    mae_test        = mae_test,
    benchmark       = benchmark_tbl,
    best_model      = best_model_info,
    stratification  = stratification,
    km_plot         = km_plot,
    heatmap         = heatmap_plot,
    c_index_heatmap = c_index_heatmap,
    enrichment      = enrichment_res
  )
}


# =============================================================================
# Example usage
# =============================================================================
if (FALSE) {

  #  1: download from LinkedOmics -------------------------
  results_lgg <- run_Omics2Surv(
    cancer_type       = "TCGA-KIPAN",
    linkedomics_omics = list(
      mRNA   = "RNASeq (HiSeq, Gene Level)",
      methyl = "Methylation (Gene level, HM450K)"
    ),
    screening_specs = list(
      mRNA   = list(type = "RNA",   method = "supervised",
                    supervised_method = "pval", p_threshold = 0.05,
                    q_threshold = 0.25),
      methyl = list(type = "methyl", method = "unsupervised",
                    unsupervised_method = "IQR", hvg_threshold = 25)
    ),
    model_specs = list(
      list(type = "pen_reg", omic_list = "mRNA",
           integration = "none", model_method = "cox_lasso"),
      list(type = "pen_reg", omic_list = c("mRNA", "methyl"),
           integration = "early", model_method = "cox_elastic_net")
    ),
    train_fraction = 0.7,
    strat_method   = "survcutpoint",
    dataset_name   = "TCGA-KIPAN"
  )

  # 2: user data, with splitting ------------------------------------
  my_rna   <- prepare_custom_data(read.csv("rna.csv",   row.names = 1), "Omics")
  my_clin  <- prepare_custom_data(read.csv("clin.csv",  row.names = 1), "Clinical")

  results_custom <- run_Omics2Surv(
    omics_list    = list(mRNA = my_rna),
    clinical_data = my_clin,
    screening_specs = list(
      mRNA = list(type = "RNA", method = "unsupervised",
                  unsupervised_method = "MAD", top_n_features = 500)
    ),
    model_specs = list(
      list(type = "pen_reg", omic_list = "mRNA",
           integration = "none", model_method = "cox_lasso")
    ),
    train_fraction = 0.75,
    dataset_name   = "MyDataset"
  )

  # 3: no split (just training, best by Harrell train) --------
  results_notrain <- run_Omics2Surv(
    omics_list    = list(mRNA = my_rna),
    clinical_data = my_clin,
    screening_specs = list(
      mRNA = list(type = "RNA", method = "supervised",
                  supervised_method = "pval", p_threshold = 0.05)
    ),
    model_specs = list(
      list(type = "pen_reg", omic_list = "mRNA",
           integration = "none", model_method = "cox_lasso")
    ),
    train_fraction = NULL     # no split
  )

  # Inspect
  results_lgg$benchmark
  results_lgg$best_model
  print(results_lgg$km_plot)
  print(results_lgg$c_index_heatmap)   # NULL if train_fraction = NULL
}
