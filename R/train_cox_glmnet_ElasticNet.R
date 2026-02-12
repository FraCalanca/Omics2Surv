#' Train a Cox Elastic Net Model with Alpha Optimization
#'
#' Fits a Cox proportional hazards model using Elastic Net regularization. The function
#' iterates through a range of \code{alpha} values to identify the optimal mixture
#' of L1 and L2 penalties based on cross-validated performance.
#'
#' @param x_training A matrix or data frame of predictor variables.
#' @param y_training A data frame containing survival outcomes with columns
#'   `overall_survival` and `status`.
#' @param nfolds Integer. Number of folds for cross-validation. Default is 5.
#' @param optimal_lambda Character. Rule for lambda selection: \code{"minimum"}
#'   (lambda.min) or \code{"1se"} (lambda.1se).
#' @param type_measure Character. Performance metric for CV (e.g., \code{"C"}
#'   for C-index or \code{"deviance"}).
#' @param alpha_range Numeric vector. The sequence of \code{alpha} values (0 to 1)
#'   to test. Default is 0.1 to 0.9.
#' @param penalty Optional numeric vector of penalty factors for individual features.
#' @param system_type Character. Operating system: \code{"unix"} or \code{"windows"}.
#' @param nlambda Integer. Number of \eqn{\lambda} values in each \code{glmnet} path.
#'   Default is 200.
#' @param seed Integer. Seed for reproducibility.
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{cv_fit}: The \code{cv.glmnet} object for the best performing \code{alpha}.
#'   \item \code{best_alpha}: The alpha value that yielded the best CV result.
#'   \item \code{selected_lambda}: The numeric \eqn{\lambda} used for the final model.
#'   \item \code{selected_features}: Predictors with non-zero coefficients at the chosen alpha/lambda.
#'   \item \code{pred}: Linear predictors for the training set.
#'   \item \code{cindex}: Training set Concordance Index.
#' }
#' @export
train_cox_glmnet_ElasticNet = function(x_training, y_training,
                                       nfolds = 5,
                                       optimal_lambda = "minimum",
                                       type_measure = "C",
                                       alpha_range = seq(0.1, 0.9, by = 0.1),
                                       penalty = NULL,
                                       system_type = "unix",
                                       nlambda=200,
                                       seed = 123) {
  set.seed(seed)

  type_measure = match.arg(type_measure)
  optimal_lambda = match.arg(optimal_lambda)
  system_type = match.arg(system_type)

  time = as.numeric(y_training$overall_survival)
  event = as.numeric(y_training$status)
  x_training = as.matrix(x_training)
  y_surv = survival::Surv(time, event)

  num_cores = max(1, parallel::detectCores() - 2)
  message("Using backend with ", num_cores, " core(s).")

  if (system_type == "unix") {
    doParallel::registerDoParallel(cores = num_cores)
  } else if (system_type == "windows") {
    cl <- parallel::makeCluster(num_cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))  # stop cluster automatically when function exits
  } else {
    stop("Error: Specify 'unix' or 'windows' in system_type.")
  }

  # Perform cross-validation for different alpha values
  cv_results = lapply(alpha_range, function(a) {
    message("Testing alpha = ", a)

    if (is.null(penalty)) {
      cv_fit = glmnet::cv.glmnet(x_training, y_surv, family = "cox",
                                 alpha = a, nfolds = nfolds,nlambda=nlambda,
                                 type.measure = type_measure, parallel = TRUE, maxit=1e+09)
    } else {
      cv_fit = glmnet::cv.glmnet(x_training, y_surv, family = "cox",
                                 alpha = a, nfolds = nfolds,nlambda=200,
                                 penalty.factor = penalty,
                                 type.measure = type_measure, parallel = TRUE, maxit=1e+09)
    }

    selected_lambda = if (optimal_lambda == "minimum") cv_fit$lambda.min else cv_fit$lambda.1se
    list(alpha = a, cv_fit = cv_fit, selected_lambda = selected_lambda, cvm = cv_fit$cvm[cv_fit$lambda == selected_lambda])
  })

  # Select best alpha (minimum cross-validation error)
  if (type_measure=="deviance"){
    best_result = cv_results[[which.min(sapply(cv_results, function(res) res$cvm))]]
  } else {
    best_result = cv_results[[which.max(sapply(cv_results, function(res) res$cvm))]]
  }
  best_alpha = best_result$alpha
  best_cv_fit = best_result$cv_fit
  selected_lambda = best_result$selected_lambda

  message("Best alpha selected: ", best_alpha)

  coefficients = if (optimal_lambda == "minimum") coef(best_cv_fit, s="lambda.min") else coef(best_cv_fit, s="lambda.1se")
  selected_features = rownames(coefficients)[coefficients[, 1] != 0]

  pred = if (optimal_lambda == "minimum") {
    predict(best_cv_fit, newx = x_training, type = "link", s = "lambda.min")
  } else{
    predict(best_cv_fit, newx = x_training, type = "link", s = "lambda.1se")
  }
  train_cindex = survival::concordance(y_surv ~ pred, reverse = TRUE)

  return(list(
    cv_fit = best_cv_fit,
    best_alpha = best_alpha,
    selected_lambda = selected_lambda,
    selected_features = selected_features,
    pred=pred,
    cindex = train_cindex$concordance
  ))
}
