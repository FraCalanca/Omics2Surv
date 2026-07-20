#' Train a Penalized Accelerated Failure Time (AFT) Model
#'
#' Fits a penalized AFT model using the \code{penAFT} package. The function
#' optimizes the \code{alpha} parameter through cross-validation to balance
#' regularization types, supporting both standard Elastic Net and Sparse Group LASSO.
#'
#' @param x_training A matrix or data frame of predictors (samples in rows).
#' @param y_training A data frame containing survival data with \code{overall_survival}
#'   and \code{status} columns.
#' @param nfolds Integer. Number of folds for cross-validation. Default is 5.
#' @param penalty Character. The penalty type to apply: \code{"EN"} (Elastic Net)
#'   or \code{"SG"} (Sparse Group LASSO). Default is \code{"EN"}.
#' @param alpha_range Numeric vector. The range of \code{alpha} values to test
#'   (mixture between L1 and L2, or between Group and individual LASSO).
#' @param groups A vector of integers indicating group membership for predictors.
#'   Required if \code{penalty = "SG"}.
#' @param nlambda Integer. Number of \eqn{\lambda} values in the tuning grid.
#' @param seed Integer or NULL. Seed for reproducibility.
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{best_cv_fit}: The \code{penAFT.cv} object for the optimal \code{alpha}.
#'   \item \code{best_alpha}: The selected \code{alpha} value.
#'   \item \code{selected_lambda}: The value of \eqn{\lambda} that minimized CV error.
#'   \item \code{selected_features}: Names of features with non-zero coefficients.
#'   \item \code{coefficients}: Vector of estimated coefficients.
#'   \item \code{pred}: Linear predictors for the training set.
#'   \item \code{cv_results}: A list of results for all tested \code{alpha} values.
#' }
#' @export
train_penAFT= function(x_training, y_training,
                       nfolds = 5,
                       penalty = c("EN", "SG"),
                       alpha_range = seq(0.1, 0.9, by = 0.1),
                       groups = NULL,
                       nlambda = 100,
                       seed = 123) {
  # ---- Setup ----
  set.seed(seed)
  penalty = match.arg(penalty)
  #input checks
  if (penalty == "SG" && is.null(groups)) {
    stop("You must specify 'groups' when penalty = 'SG'.")
  }
  # ---- data preprocessing ----
  time = as.numeric(y_training$overall_survival)
  event = as.numeric(y_training$status)
  x_training = as.matrix(x_training)
  y_surv = survival::Surv(time, event)
  # ---- alpha ----
  message("Starting alpha search...")
  cv_results = lapply(alpha_range, function(a) {
    message("Testing alpha = ", a)
    aft_cv = penAFT::penAFT.cv(
      X = x_training,
      logY = log(time),
      delta = event,
      nfolds = nfolds,
      nlambda = nlambda,
      penalty = penalty,
      alpha = a,
      groups = groups,
      standardize = TRUE,
      quiet = TRUE
    )
    # get predictions for training data
    pred = as.vector(penAFT::penAFT.predict(aft_cv, x_training))
    cindex = survival::concordance(y_surv ~ pred)
    list(alpha = a, aft_cv = aft_cv, cindex = cindex)
  })
  # remove NULL (if any)
  cv_results = Filter(Negate(is.null), cv_results)
  # ---- select best alpha ----
  best_result = cv_results[[which.max(sapply(cv_results, function(res) res$cindex$concordance))]]
  best_alpha = best_result$alpha
  best_cv_fit = best_result$aft_cv
  message("Best alpha selected: ", best_alpha)
  # ---- coefficients ----
  coef_final = penAFT::penAFT.coef(best_cv_fit)
  coef_final = as.vector(coef_final$beta)
  names(coef_final) = colnames(x_training)
  selected_features = names(coef_final[coef_final != 0])
  # ---- predictions ----
  linear_pred = as.vector(penAFT::penAFT.predict(best_cv_fit, x_training))
  names(linear_pred) <- rownames(x_training)
  # ---- Output ----
  return(list(
    best_cv_fit = best_cv_fit,
    best_alpha = best_alpha,
    selected_lambda = best_cv_fit$lambda.min,
    selected_features = selected_features,
    coefficients = coef_final,
    pred = linear_pred,
    cv_results = cv_results
  ))
}


#' Evaluate a Penalized AFT Model on Test Data
#'
#' Computes linear predictors for a fitted penalized AFT model on independent
#' test data, using the cross-validated fit produced by \code{\link{train_penAFT}}.
#'
#' @param x_test A matrix or data frame of predictors for the test set
#'   (samples in rows, same features as the training set).
#' @param training_result A list returned by \code{\link{train_penAFT}}. Must
#'   contain \code{best_cv_fit}, the \code{penAFT.cv} object used for prediction.
#'
#' @return A named list with:
#' \itemize{
#'   \item \code{pred}: Linear predictors for the test set, named by the row
#'     names of \code{x_test}.
#' }
#'
#' @seealso \code{\link{train_penAFT}}
#'
#' @export
test_penAFT = function(x_test, training_result) {

  X_test <- as.matrix(x_test)

  # --- penAFT.predict ---
  linear_pred <- as.vector(
    penAFT::penAFT.predict(
      training_result$best_cv_fit,
      Xnew = X_test
    )
  )
  names(linear_pred) <- rownames(x_test)

  # --- Output ---
  return(list(
    pred = linear_pred
  ))
}
