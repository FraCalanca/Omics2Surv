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
#' @param type_measure Character. Loss function for cross-validation:
#'   \code{"C"} (C-index, default) or \code{"deviance"}.
#' @param penalty Optional numeric vector for penalty factors.
#' @param system_type Character. Operating system type, either \code{"unix"}
#'   or \code{"windows"}, to handle parallel backend registration.
#' @param nlambda Integer. The number of \eqn{\lambda} values in the sequence.
#'   Default is 100.
#' @param seed Integer. Seed for reproducibility of cross-validation folds.
#' @param ... Additional arguments (currently unused).
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{optimal_lambda}: The selection rule used.
#'   \item \code{type_measure}: The cross-validation loss function used.
#'   \item \code{cv_fit}: The fitted \code{cv.glmnet} object.
#'   \item \code{selected_lambda}: The numeric value of the chosen \eqn{\lambda}.
#'   \item \code{selected_features}: Character vector of predictors with non-zero coefficients.
#'   \item \code{pred}: Linear predictors for the training set.
#' }
#' @export
train_cox_glmnet_LASSO = function(x_training, y_training,
                                  nfolds = 5,
                                  optimal_lambda = c("minimum", "1se"),
                                  type_measure = c("C", "deviance"),
                                  penalty = NULL,
                                  system_type = "unix",
                                  nlambda=100,
                                  seed = 123,
                                  ...) {

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
                               type.measure = type_measure, parallel = TRUE, maxit=1e+05, cox.ties = "efron")
  } else {
    cv_fit = glmnet::cv.glmnet(x_training, y_surv, family = "cox",
                               alpha = 1, nfolds = nfolds,
                               penalty.factor = penalty,
                               type.measure = type_measure, parallel = TRUE, maxit=1e+06, cox.ties = "efron")
  }

  selected_lambda = if (optimal_lambda == "minimum") cv_fit$lambda.min else cv_fit$lambda.1se

  coefficients = if (optimal_lambda == "minimum") coef(cv_fit, s="lambda.min") else coef(cv_fit, s="lambda.1se")
  selected_features = rownames(coefficients)[coefficients[, 1] != 0]

  pred = if (optimal_lambda == "minimum") {
    predict(cv_fit, newx = x_training, type = "link", s = "lambda.min")
    } else{
      predict(cv_fit, newx = x_training, type = "link", s = "lambda.1se")
    }

  return(list(optimal_lambda = optimal_lambda,
              type_measure = type_measure,
              cv_fit = cv_fit,
              selected_lambda = selected_lambda,
              selected_features = selected_features,
              pred = pred))
}


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
#' @param type_measure Character. Performance metric for CV: \code{"C"}
#'   (C-index, default) or \code{"deviance"}.
#' @param alpha_range Numeric vector. The sequence of \code{alpha} values (0 to 1)
#'   to test. Default is 0.1 to 0.9.
#' @param penalty Optional numeric vector of penalty factors for individual features.
#' @param system_type Character. Operating system: \code{"unix"} or \code{"windows"}.
#' @param nlambda Integer. Number of \eqn{\lambda} values in each \code{glmnet} path.
#'   Default is 100.
#' @param seed Integer. Seed for reproducibility.
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{cv_fit}: The \code{cv.glmnet} object for the best performing \code{alpha}.
#'   \item \code{best_alpha}: The alpha value that yielded the best CV result.
#'   \item \code{selected_lambda}: The numeric \eqn{\lambda} used for the final model.
#'   \item \code{selected_features}: Predictors with non-zero coefficients at the chosen alpha/lambda.
#'   \item \code{pred}: Linear predictors for the training set.
#' }
#' @export
train_cox_glmnet_ElasticNet = function(x_training, y_training,
                                       nfolds = 5,
                                       optimal_lambda = c("minimum", "1se"),
                                       type_measure = c("C", "deviance"),
                                       alpha_range = seq(0.1, 0.9, by = 0.1),
                                       penalty = NULL,
                                       system_type = c("unix", "windows"),
                                       nlambda=100,
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
                                 type.measure = type_measure, parallel = TRUE, maxit=1e+06, cox.ties = "efron")
    } else {
      cv_fit = glmnet::cv.glmnet(x_training, y_surv, family = "cox",
                                 alpha = a, nfolds = nfolds,nlambda=100,
                                 penalty.factor = penalty,
                                 type.measure = type_measure, parallel = TRUE, maxit=1e+06,cox.ties = "efron")
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

  return(list(
    cv_fit = best_cv_fit,
    best_alpha = best_alpha,
    selected_lambda = selected_lambda,
    selected_features = selected_features,
    pred=pred
  ))
}


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
#' @param type_measure Character. Loss function for cross-validation:
#'   \code{"C"} (C-index, default) or \code{"deviance"}.
#' @param unpenalized Optional vector of indices or names for variables to be
#'   excluded from the penalty (penalty set to 0).
#' @param system_type Character. OS type: \code{"unix"} or \code{"windows"}.
#' @param nlambda Integer. Number of \eqn{\lambda} values for the Ridge step sequence.
#' @param gamma Numeric. Relaxation parameter passed to \code{glmnet} for the
#'   Ridge step. Default is 1.
#' @param seed Optional integer for reproducibility.
#' @param ... Additional arguments (currently unused).
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{ridge_cv}: The preliminary Ridge cross-validation object.
#'   \item \code{cv_fit}: The final Adaptive LASSO cross-validation object.
#'   \item \code{selected_lambda}: The numeric \eqn{\lambda} used for the final model.
#'   \item \code{selected_features}: Predictors with non-zero coefficients.
#'   \item \code{pred}: Linear predictors for the training set.
#'   \item \code{weights}: The penalty factors derived from the Ridge step.
#' }
#' @export
train_cox_glmnet_adaptiveLASSO = function(x_training, y_training,
                                          nfolds = 5,
                                          optimal_lambda = c("minimum", "1se"),
                                          type_measure = c("C", "deviance"),
                                          unpenalized = NULL,
                                          system_type = "unix",
                                          nlambda=100,
                                          gamma = 1,
                                          seed = NULL,
                                          ...) {
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
                               gamma = gamma,
                               parallel = TRUE, maxit = 1e+06, nlambda= nlambda, trace.it=1, cox.ties = "efron")

  ridge_coef = as.vector(coef(ridge_cv, s = "lambda.1se"))
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
                             nlambda = nlambda,
                             penalty.factor = weights,
                             type.measure = type_measure,
                             parallel = TRUE, maxit = 1e+06, trace.it=1, cox.ties = "efron")

  selected_lambda = if (optimal_lambda == "minimum") cv_fit$lambda.min else cv_fit$lambda.1se

  coefficients = if (optimal_lambda == "minimum") coef(cv_fit, s="lambda.min") else coef(cv_fit, s="lambda.1se")
  selected_features = rownames(coefficients)[coefficients[, 1] != 0]

  pred = if (optimal_lambda == "minimum") {
    predict(cv_fit, newx = x_training, type = "link", s = "lambda.min")
  } else{
    predict(cv_fit, newx = x_training, type = "link", s = "lambda.1se")
  }

  return(list(ridge_cv = ridge_cv,
              cv_fit = cv_fit,
              selected_lambda = selected_lambda,
              selected_features = selected_features,
              pred = pred,
              weights = weights))
}


#' Evaluate Cox GLMNET Model on Test Data
#'
#' Predicts risk scores for a fitted Cox model using independent test data.
#'
#' @param x_test A matrix or data frame of predictor variables for the test set.
#' @param training_result A list containing the trained model objects:
#'   \itemize{
#'     \item \code{cv_fit}: A fitted \code{cv.glmnet} object.
#'     \item \code{selected_lambda}: The penalty parameter (\eqn{\lambda}) to use for prediction.
#'   }
#'
#' @return A list with predictions:
#' \code{pred}: A matrix of predicted risk scores (linear predictors) for the test set.
#'
#' @seealso \code{\link{train_cox_glmnet_LASSO}},
#'   \code{\link{train_cox_glmnet_ElasticNet}},
#'   \code{\link{train_cox_glmnet_adaptiveLASSO}}
#'
#' @export
test_cox_glmnet= function (x_test, training_result) {

  x_test = as.matrix(x_test)

  pred_test= predict (training_result$cv_fit, newx = x_test, type = "link", s = training_result$selected_lambda)

  return(list(pred = pred_test))
}
