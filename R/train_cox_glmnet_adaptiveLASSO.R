#' Train an Adaptive Cox LASSO Model
#'
#' Fits a Cox proportional hazards model using Adaptive LASSO. This two-step
#' procedure uses Ridge regression coefficients to create penalty factors,
#' encouraging the model to select variables with stronger initial signals.
#'
#' @param x_training A matrix or data frame of predictor variables.
#' @param y_training A data frame containing survival outcomes, specifically
#'   `overall_survival` (time) and `status` (event).
#' @param nfolds Integer. Number of folds for cross-validation in both steps.
#' @param optimal_lambda Character. Choice of lambda for prediction:
#'   \code{"minimum"} or \code{"1se"}.
#' @param type_measure Character. Loss function for cross-validation (e.g., \code{"C"}).
#' @param unpenalized Optional vector of indices or names for variables to be
#'   excluded from the penalty (penalty set to 0).
#' @param system_type Character. OS type: \code{"unix"} or \code{"windows"}.
#' @param nlambda Integer. Number of \eqn{\lambda} values for the Ridge step sequence.
#' @param seed Optional integer for reproducibility.
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{ridge_cv}: The preliminary Ridge cross-validation object.
#'   \item \code{cv_fit}: The final Adaptive LASSO cross-validation object.
#'   \item \code{selected_lambda}: The numeric \eqn{\lambda} used for the final model.
#'   \item \code{selected_features}: Predictors with non-zero coefficients.
#'   \item \code{pred}: Linear predictors for the training set.
#'   \item \code{cindex}: Training set Concordance Index.
#'   \item \code{weights}: The penalty factors derived from the Ridge step.
#' }
#' @export
train_cox_glmnet_adaptiveLASSO = function(x_training, y_training,
                                          nfolds = 5,
                                          optimal_lambda = "minimum",
                                          type_measure = "C",
                                          unpenalized = NULL,
                                          system_type = "unix",
                                          nlambda=150,
                                          seed = NULL) {
  if (!is.null(seed)) {set.seed(seed)}

  type_measure = match.arg(type_measure)
  optimal_lambda = match.arg(optimal_lambda)

  time = as.numeric(y_training$overall_survival)
  event = as.numeric(y_training$status)
  x_training = as.matrix(x_training)
  y_surv = survival::Surv(time, event)

  num_cores = max(1, parallel::detectCores() -2)
  message("Using ", num_cores, " core(s).")

  if (system_type == "unix") {
    doParallel::registerDoParallel(cores = num_cores)
  } else if (system_type == "windows") {
    cl <- parallel::makeCluster(num_cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  } else {
    stop("Specify 'unix' or 'windows'.")
  }

  # --- Step 1: Ridge Cox for weights ---
  ridge_cv = glmnet::cv.glmnet(x_training, y_surv,
                               family = "cox",
                               alpha = 0,
                               nfolds = nfolds,
                               type.measure = type_measure,
                               parallel = TRUE, maxit = 1e+09, nlambda= nlambda, trace.it=1)

  ridge_coef = as.vector(coef(ridge_cv, s = "lambda.min"))
  weights = abs(1 / ridge_coef)

  # handle infinite weights (e.g. zero coefficients)
  weights[!is.finite(weights)] = max(weights[is.finite(weights)], na.rm = TRUE)

  # --- Optionally exclude some variables from penalization ---
  if (!is.null(unpenalized)) {
    weights[unpenalized] = 0
  }

  # --- Step 2: Adaptive LASSO Cox ---
  cv_fit = glmnet::cv.glmnet(x_training, y_surv,
                             family = "cox",
                             alpha = 1,
                             nfolds = nfolds,
                             penalty.factor = weights,
                             type.measure = type_measure,
                             parallel = TRUE, maxit = 1e+09, trace.it=1)

  selected_lambda = if (optimal_lambda == "minimum") cv_fit$lambda.min else cv_fit$lambda.1se

  coefficients = if (optimal_lambda == "minimum") coef(cv_fit, s="lambda.min") else coef(cv_fit, s="lambda.1se")
  selected_features = rownames(coefficients)[coefficients[, 1] != 0]

  pred = if (optimal_lambda == "minimum") {
    predict(cv_fit, newx = x_training, type = "link", s = "lambda.min")
  } else{
    predict(cv_fit, newx = x_training, type = "link", s = "lambda.1se")
  }
  train_cindex = survival::concordance(y_surv ~ pred, reverse = TRUE)

  return(list(ridge_cv = ridge_cv,
              cv_fit = cv_fit,
              selected_lambda = selected_lambda,
              selected_features = selected_features,
              pred = pred,
              cindex = train_cindex$concordance,
              weights = weights))
}
