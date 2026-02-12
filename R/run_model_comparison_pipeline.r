#' Run Model Comparison Pipeline for Multi-Omics Survival Analysis
#'
#' Executes a comprehensive comparison of multiple survival models across different
#' integration strategies (none, early, late, or joint), returning performance metrics
#' for each model.
#'
#' @param omics_list Named list of SummarizedExperiment objects, each representing
#'   an omics layer.
#' @param clinical_data Matrix or data frame of clinical data (variables as rows,
#'   samples as columns).
#' @param surv_list Character vector of survival variable names (e.g.,
#'   \code{c("overall_survival", "status")}).
#' @param integration Character string specifying integration strategy: \code{"none"}
#'   (single omics), \code{"early"} (concatenated features), \code{"late"} (ensemble
#'   predictions), or \code{"joint"} (cooperative modeling).
#' @param model_methods Character vector of model names to compare. Options depend on
#'   \code{integration}: for "none"/"early", use penalized regression models; for
#'   "late", use same; for "joint", use \code{c("CoopCox", "AFTCoop", "blockForest",
#'   "flexynesis")}.
#' @param train_fraction Numeric proportion (0-1) of samples for training. Default is 0.7.
#' @param seed Integer seed for reproducibility. Default is 123.
#' @param late_method Character string for late integration: \code{"average"} or
#'   \code{"RRA"}. Only used when \code{integration = "late"}. Default is \code{"average"}.
#' @param ... Additional arguments passed to model-specific functions.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{comparison}: Data frame with columns \code{model}, \code{train_cindex},
#'       and \code{test_cindex} for each model tested
#'     \item \code{per_model_results}: Named list of detailed results for each model,
#'       including reports from \code{metadata(mae)$analysis}
#'   }
#'
#' @export
run_model_comparison_pipeline <- function(
    omics_list,
    clinical_data,
    surv_list,
    integration = c("none","early","late","joint"),
    model_methods,
    train_fraction = 0.7,
    seed = 123,
    late_method = c("average","RRA"),
    ...
) {

  integration <- match.arg(integration)
  late_method <- match.arg(late_method)

  # -------------------------
  # 0. Filter common patients
  # -------------------------
  omics_list <- filter_common_patients_mae(omics_list)

  # -------------------------
  # 1. Build MAE
  # -------------------------
  mae <- build_MAE(omics_list, clinical_data)

  # -------------------------
  # 2. Split
  # -------------------------
  assay_for_split <- names(omics_list)[1]

  X <- prepare_design_matrix(mae, assay_name = assay_for_split)
  y <- prepare_surv_data(mae, surv_list, sample_order = rownames(X))

  split <- split_train_test(
    X, y,
    train_fraction = train_fraction,
    seed = seed
  )

  # -------------------------
  # 3. Loop on models
  # -------------------------
  results <- list()

  for (m in model_methods) {

    mae_tmp <- mae

    if (integration %in% c("none","early")) {

      mae_tmp <- run_survival_model_on_mae(
        mae_tmp,
        assay_name = assay_for_split,
        split = split,
        method = m,
        seed = seed,
        ...
      )

      rep <- metadata(mae_tmp)$analysis[[m]][[assay_for_split]]$report
      results[[m]] <- rep
    }

    if (integration == "late") {

      mae_tmp <- run_late_integration_on_mae(
        mae_tmp,
        omic_list = names(omics_list),
        surv_list = surv_list,
        model_method = m,
        train_fraction = train_fraction,
        seed = seed,
        late_method = late_method,
        ...
      )

      rep <- metadata(mae_tmp)$analysis$late_integration[[m]]$report
      results[[m]] <- rep
    }

    if (integration == "joint") {

      mae_tmp <- run_joint(
        mae = mae_tmp,
        omic_list = names(omics_list),
        surv_list = surv_list,
        joint_method = m,
        train_fraction = train_fraction,
        seed = seed,
        ...
      )

      rep <- metadata(mae_tmp)$analysis$joint[[m]]$report
      results[[m]] <- rep
    }
  }

  # -------------------------
  # 4. Comparison table
  # -------------------------
  comparison_table <- do.call(
    rbind,
    lapply(names(results), function(m) {
      data.frame(
        model = m,
        train_cindex = results[[m]]$train_cindex %||% NA,
        test_cindex  = results[[m]]$test_cindex  %||% results[[m]]$consensus_cindex
      )
    })
  )

  rownames(comparison_table) <- NULL

  list(
    comparison = comparison_table,
    per_model_results = results
  )
}
