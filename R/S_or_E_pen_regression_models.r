#' Single or Early Integration Penalized Regression Models
#'
#' Trains and tests penalized regression survival models on a single omics layer or
#' early-integrated (concatenated) features, storing results in the MultiAssayExperiment
#' metadata.
#'
#' @param mae A MultiAssayExperiment object.
#' @param assay_name Character string specifying the assay name to model (e.g., single
#'   omics layer or "early_integrated").
#' @param split A list containing:
#'   \itemize{
#'     \item \code{train_omic}: Training design matrix (samples × features)
#'     \item \code{test_omic}: Test design matrix
#'     \item \code{train_survival}: Training survival data frame with
#'       \code{overall_survival} and \code{status}
#'     \item \code{test_survival}: Test survival data frame
#'   }
#' @param method Character string specifying the penalized model: \code{"cox_lasso"},
#'   \code{"cox_adaptive_lasso"}, \code{"cox_elastic_net"}, or \code{"aft_pen"}.
#' @param seed Integer seed for reproducibility. Default is 123.
#' @param ... Additional arguments passed to the training function (e.g., \code{alpha},
#'   \code{nfolds}).
#'
#' @return The input MultiAssayExperiment object with results stored at
#'   \code{metadata(mae)$analysis[[method]][[assay_name]]} containing:
#'   \itemize{
#'     \item \code{train}: Training results (model, predictions, C-index, selected features)
#'     \item \code{test}: Test results (predictions, C-index)
#'   }
#'
#' @export
S_or_E_pen_regression_models <- function(
    mae,
    assay_name,
    split,
    method = c("cox_lasso","cox_adaptive_lasso","cox_elastic_net",
               "aft_pen"),
    seed = 123,
    ...
){

  method <- match.arg(method)

  x_train <- split$train_omic
  y_train <- split$train_survival
  x_test  <- split$test_omic
  y_test  <- split$test_survival

  # train/test
  train_fun <- switch(method,
                      cox_lasso = train_cox_glmnet_LASSO,
                      cox_adaptive_lasso = train_cox_glmnet_adaptiveLASSO,
                      cox_elastic_net = train_cox_glmnet_ElasticNet,
                      aft_pen = train_penAFT)

  test_fun <- switch(method,
                     cox_lasso = test_cox_glmnet,
                     cox_adaptive_lasso = test_cox_glmnet,
                     cox_elastic_net = test_cox_glmnet,
                     aft_pen = test_penAFT)

  # -------------------------
  # 1. Training
  # -------------------------
  train_res <- train_fun(x_train, y_train, seed = seed, ...)

  # -------------------------
  # 2. Test
  # -------------------------
  test_res <- test_fun(x_test, y_test, train_res)

  # -------------------------
  # 3. Recording in MAE
  # -------------------------
  if (is.null(metadata(mae)$analysis)) metadata(mae)$analysis <- list()
  if (is.null(metadata(mae)$analysis[[method]])) metadata(mae)$analysis[[method]] <- list()

  metadata(mae)$analysis[[method]][[assay_name]] <- list(
    train = train_res,
    test  = test_res
  )

  mae
}
