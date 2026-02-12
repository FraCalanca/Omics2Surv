#' Evaluate Cox GLMNET Model on Test Data
#'
#' Predicts risk scores and calculates the Concordance Index (C-index) for a
#' fitted Cox lasso model using independent test data.
#'
#' @param x_test A matrix or data frame of predictor variables for the test set.
#' @param y_test A data frame or matrix containing survival information.
#'   Must include columns `overall_survival` (time) and `status` (event).
#' @param training_result A list containing the trained model objects:
#'   \itemize{
#'     \item \code{cv_fit}: A fitted \code{cv.glmnet} object.
#'     \item \code{selected_lambda}: The penalty parameter (\eqn{\lambda}) to use for prediction.
#'   }
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{pred}: A matrix of predicted risk scores (linear predictors) for the test set.
#'   \item \code{cindex}: A numeric value representing the Harrell's Concordance Index.
#' }
#' @export
test_cox_glmnet= function (x_test, y_test, training_result) {

  x_test = as.matrix(x_test)
  time = as.numeric(y_test$overall_survival)
  event = as.numeric(y_test$status)
  y_test_surv = survival::Surv(time, event)

  pred_test= predict (training_result$cv_fit, newx = x_test, type = "link", s = training_result$selected_lambda)
  Cindex_test= survival::concordance(y_test_surv ~ pred_test, reverse = TRUE)

  return(list(pred = pred_test, cindex = Cindex_test$concordance))
}
