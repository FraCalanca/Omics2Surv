#' Run Penalized Regression Models with Multi-Omics Integration
#'
#' Applies penalized regression survival models (Lasso, adaptive Lasso, elastic net,
#' or AFT) using single omics, early integration (concatenated features), or late
#' integration (ensemble predictions).
#'
#' @param omics_list Named list of SummarizedExperiment objects, each representing
#'   an omics layer.
#' @param clinical_data Matrix or data frame of clinical data (variables as rows,
#'   samples as columns).
#' @param surv_list Character vector of survival variable names (e.g.,
#'   \code{c("overall_survival", "status")}).
#' @param integration Character string specifying integration strategy: \code{"none"}
#'   (single omics layer), \code{"early"} (concatenate features before modeling),
#'   or \code{"late"} (model each layer separately, then ensemble).
#' @param model_method Character string specifying the penalized model: \code{"cox_lasso"},
#'   \code{"cox_adaptive_lasso"}, \code{"cox_elastic_net"}, or \code{"aft_pen"}.
#' @param train_fraction Numeric proportion (0-1) of samples for training. Default is 0.7.
#' @param seed Integer seed for reproducibility. Default is 123.
#' @param late_method Character string for late integration: \code{"average"} or
#'   \code{"RRA"}. Only used when \code{integration = "late"}. Default is \code{"average"}.
#' @param ... Additional arguments passed to the penalized regression functions.
#'
#' @return A MultiAssayExperiment object with results stored in \code{metadata(mae)$analysis}:
#'   \itemize{
#'     \item For \code{"none"} or \code{"early"}: Results at
#'       \code{metadata(mae)$analysis[[model_method]][[assay_name]]}
#'     \item For \code{"late"}: Results at
#'       \code{metadata(mae)$analysis$late_integration[[model_method]]}
#'   }
#'   Each contains \code{train} and \code{test} results with predictions and C-index.
#'
#' @export
run_penalized_regression_models <- function(
    omics_list,
    clinical_data,
    surv_list,
    integration = c("none", "early", "late"),
    model_method = c("cox_lasso","cox_adaptive_lasso","cox_elastic_net",
                     "aft_pen"),
    train_fraction = 0.7,
    seed = 123,
    late_method = c("average","RRA"),
    ...
) {

  integration <- match.arg(integration)
  model_method <- match.arg(model_method)
  late_method <- match.arg(late_method)

  # 1. MAE
  Mae <- build_MAE(omics_list, clinical_data)

  # 2. Single / early integration
  if (integration == "early") {
    omics_list<- filter_common_patients(omics_list)
    Mae <- early_integration(Mae, names(omics_list))
    assay_name <- "early_integrated"

  } else if (integration == "none") {
    assay_name <- names(omics_list)
  }

  if (integration %in% c("none","early")) {
    X <- prepare_design_matrix(Mae, assay_name = assay_name)
    y <- prepare_surv_data(Mae, surv_list, sample_order = rownames(X))
    split <- split_train_test(X, y, train_fraction = train_fraction, seed = seed)

    Mae <- S_or_E_pen_regression_models(Mae, assay_name = assay_name, split = split,
                                     method = model_method, seed = seed, ...)
  }

  # 4. Late integration
  if (integration == "late") {
    Mae <- L_pen_regression_models(mae = Mae, omic_list = omics_list, surv_list = surv_list,
                                   method = late_method, model = model_method, train_fraction = train_fraction,
                                   seed = seed, ... )
  }

  Mae
}
