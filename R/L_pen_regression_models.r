#' Late Integration of Penalized Regression Models Across Multiple Omics Layers
#'
#' Performs late integration of survival predictions from penalized regression
#' models trained independently on each omics layer. Each omics dataset is modeled
#' separately using the same train/test split, and predictions are combined using
#' averaging or rank-based aggregation methods.
#'
#' @param mae A MultiAssayExperiment object containing multiple omics datasets and
#'   clinical information. Results will be stored in
#'   \code{metadata(mae)$analysis$late_integration}.
#' @param omic_list A list of omics layer names (or SummarizedExperiment objects)
#'   corresponding to assays in the MAE. Each will be modeled independently. Common
#'   samples across all layers are automatically identified and used.
#' @param surv_list A specification for extracting survival data from the MAE,
#'   typically containing column names for time and event variables. Passed to
#'   \code{prepare_surv_data}.
#' @param method Character string specifying the late integration method for combining
#'   predictions across omics layers. Must be one of:
#'   \itemize{
#'     \item \code{"average"}: Simple or weighted average of risk scores
#'     \item \code{"RRA"}: Robust Rank Aggregation, combines rankings rather than raw scores (default)
#'   }
#'   Default is \code{"average"}.
#' @param model Character string specifying the penalized regression model to fit
#'   on each omics layer. Must be one of:
#'   \itemize{
#'     \item \code{"cox_lasso"}: Cox regression with Lasso penalty
#'     \item \code{"cox_adaptive_lasso"}: Cox regression with adaptive Lasso penalty
#'     \item \code{"cox_elastic_net"}: Cox regression with elastic net penalty
#'     \item \code{"aft_pen"}: Accelerated Failure Time model with elastic net or grouped lasso penalties
#'   }
#'   Default is \code{"cox_lasso"}.
#' @param train_fraction Numeric value between 0 and 1 specifying the proportion of
#'   samples to use for training. The remainder is used for testing. Default is 0.7
#'   (70\% train, 30\% test).
#' @param seed Integer seed for reproducibility of the train/test split. Default is 123.
#'   The same seed ensures identical splits across all omics layers.
#'
#' @return The input MultiAssayExperiment object with updated metadata. Results are
#'   stored at \code{metadata(mae)$analysis$late_integration[[model]]} containing:
#'   \itemize{
#'     \item \code{train}: Late integration results on training data, including:
#'       \itemize{
#'         \item Ensemble predictions (combined risk scores)
#'         \item Ensemble C-index
#'         \item Individual omics-specific predictions and C-indices
#'         \item Integration method used
#'       }
#'     \item \code{test}: Late integration results on test data with the same structure
#'   }
#'   Additionally, individual omics results are stored at:
#'   \code{metadata(mae)$analysis[[model]][[omic_name]]}
#'
#' @details
#' The function implements a late integration (late fusion) workflow:
#' \enumerate{
#'   \item Identifies common samples across all omics layers
#'   \item Creates a single reference train/test split based on the first omics layer
#'   \item Applies the same split to all other omics layers (ensures consistency)
#'   \item Trains the specified penalized model independently on each omics layer
#'   \item Combines predictions from all layers using the specified integration method
#'   \item Evaluates ensemble performance on both training and test sets
#'   \item Stores individual and ensemble results in MAE metadata
#' }
#'
#' @section Integration Methods:
#' \describe{
#'   \item{average}{Combines risk scores by averaging. Can be simple mean or weighted
#'     based on layer-specific performance. Assumes risk scores are on comparable scales.}
#'   \item{RRA}{Robust Rank Aggregation converts predictions to ranks within each layer,
#'     then combines ranks. More robust to outliers and scale differences. Useful when
#'     absolute risk scores vary widely across layers.}
#' }
#'
#' @section Result Organization:
#' Results are stored hierarchically:
#' \preformatted{
#' metadata(mae)
#'   └─ analysis
#'      ├─ [model]                     # Individual omics results
#'      │  ├─ [omic_1]
#'      │  │  ├─ train
#'      │  │  └─ test
#'      │  └─ [omic_2]
#'      │     ├─ train
#'      │     └─ test
#'      └─ late_integration            # Ensemble results
#'         └─ [model]
#'            ├─ train (ensemble)
#'            └─ test (ensemble)
#' }
#'
#' @seealso
#' \code{\link{S_or_E_pen_regression_models}} for single omics penalized regression
#' \code{\link{late_integration}} for the integration function
#' \code{\link{split_train_test}} for data splitting
#' \code{\link{prepare_design_matrix}} for omics data preparation
#'
#' @export
L_pen_regression_models <- function(
    mae,
    omic_list,
    surv_list,
    method = "RRA",
    model = c("cox_lasso","cox_adaptive_lasso","cox_elastic_net","aft_pen"),
    train_fraction = 0.7,
    seed = 123
) {

  method <- match.arg(method)
  model  <- match.arg(model)

  # --------------------------------------------------
  # 0. Reference split (ONCE)
  # --------------------------------------------------
  omic_list<- filter_common_patients(omic_list)

  ref_omic <- names(omic_list)[1]

  X_ref <- prepare_design_matrix(mae, assay_name = ref_omic)
  y_ref <- prepare_surv_data(mae, surv_list, sample_order = rownames(X_ref))

  split_ref <- split_train_test(
    X_ref, y_ref,
    train_fraction = train_fraction,
    seed = seed
  )

  # Containers
  res_train_list <- list()
  res_test_list  <- list()

  # --------------------------------------------------
  # 1. Loop on omics (same split)
  # --------------------------------------------------
  for (omic_name in names(omic_list)) {

    X <- prepare_design_matrix(mae, assay_name = omic_name)

    split <- list(
      train_omic     = X[rownames(split_ref$train_survival), , drop = FALSE],
      test_omic      = X[rownames(split_ref$test_survival),  , drop = FALSE],
      train_survival = split_ref$train_survival,
      test_survival  = split_ref$test_survival
    )

    mae <- S_or_E_pen_regression_models(
      mae,
      assay_name = omic_name,
      split = split,
      method = model,
      seed = seed
    )

    res_train_list[[omic_name]] <-
      metadata(mae)$analysis[[model]][[omic_name]]$train

    res_test_list[[omic_name]] <-
      metadata(mae)$analysis[[model]][[omic_name]]$test
  }

  # --------------------------------------------------
  # 2. Late integration
  # --------------------------------------------------
  late_train <- late_integration(
    res_list = res_train_list,
    y = split_ref$train_survival,
    model = model,
    method = method
  )

  late_test <- late_integration(
    res_list = res_test_list,
    y = split_ref$test_survival,
    model = model,
    method = method
  )

  # --------------------------------------------------
  # 3. Store in MAE
  # --------------------------------------------------
  if (is.null(metadata(mae)$analysis$late_integration))
    metadata(mae)$analysis$late_integration <- list()

  metadata(mae)$analysis$late_integration[[model]] <- list(
    train  = late_train,
    test   = late_test
  )

  mae
}
