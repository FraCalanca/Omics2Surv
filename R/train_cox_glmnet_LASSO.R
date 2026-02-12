#' Train a Cox LASSO Model with Cross-Validation
#'
#' Fits a Cox proportional hazards model using LASSO regularization (L1 penalty)
#' with k-fold cross-validation and parallel processing support.
#'
#' @param x_training A matrix or data frame of predictors (e.g., omics data).
#' @param y_training A data frame containing survival outcomes, specifically
#'   `overall_survival` (time) and `status` (event).
#' @param nfolds Integer. Number of folds for cross-validation. Default is 5.
#' @param optimal_lambda Character. Choice of lambda for prediction:
#'   \code{"minimum"} (lambda.min) or \code{"1se"} (lambda.1se).
#' @param type_measure Character. Loss function for cross-validation.
#'   Default is \code{"C"} (C-index).
#' @param penalty Optional numeric vector for penalty factors (e.g., for Adaptive LASSO).
#' @param system_type Character. Operating system type, either \code{"unix"}
#'   or \code{"windows"}, to handle parallel backend registration.
#' @param nlambda Integer. The number of \eqn{\lambda} values in the sequence.
#'   Default is 200.
#' @param seed Integer. Seed for reproducibility of cross-validation folds.
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{optimal_lambda}: The selection rule used.
#'   \item \code{cv_fit}: The fitted \code{cv.glmnet} object.
#'   \item \code{selected_lambda}: The numeric value of the chosen \eqn{\lambda}.
#'   \item \code{selected_features}: Character vector of predictors with non-zero coefficients.
#'   \item \code{pred}: Linear predictors for the training set.
#'   \item \code{cindex}: Training set Concordance Index.
#' }
#' @export
train_cox_glmnet_LASSO = function(x_training, y_training,
                                  nfolds = 5,
                                  optimal_lambda = "minimum",
                                  type_measure = "C",
                                  penalty = NULL,
                                  system_type = "unix",
                                  nlambda=200,
                                  seed = 123) {

  set.seed(seed)

  type_measure = match.arg(type_measure)
  optimal_lambda = match.arg(optimal_lambda)

  time = as.numeric(y_training$overall_survival)
  event = as.numeric(y_training$status)
  x_training = as.matrix(x_training)
  y_surv = survival::Surv(time, event)

  num_cores = max(1, parallel::detectCores() - 2)
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

  if (is.null(penalty)) {
    cv_fit = glmnet::cv.glmnet(x_training, y_surv, family = "cox",
                               alpha = 1, nfolds = nfolds, nlambda=nlambda,
                               type.measure = type_measure, parallel = TRUE, maxit=1e+05)
  } else {
    cv_fit = glmnet::cv.glmnet(x_training, y_surv, family = "cox",
                               alpha = 1, nfolds = nfolds,
                               penalty.factor = penalty,
                               type.measure = type_measure, parallel = TRUE, maxit=1e+09)
  }

  selected_lambda = if (optimal_lambda == "minimum") cv_fit$lambda.min else cv_fit$lambda.1se

  coefficients = if (optimal_lambda == "minimum") coef(cv_fit, s="lambda.min") else coef(cv_fit, s="lambda.1se")
  selected_features = rownames(coefficients)[coefficients[, 1] != 0]

  pred = if (optimal_lambda == "minimum") {
    predict(cv_fit, newx = x_training, type = "link", s = "lambda.min")
    } else{
      predict(cv_fit, newx = x_training, type = "link", s = "lambda.1se")
    }
  train_cindex = survival::concordance(y_surv ~ pred, reverse = TRUE)

  return(list(optimal_lambda = optimal_lambda,
              type_measure = type_measure,
              cv_fit = cv_fit,
              selected_lambda = selected_lambda,
              selected_features = selected_features,
              pred = pred,
              cindex = train_cindex$concordance))
}
