#' Train AFT Cooperative Learning Model
#'
#' Trains an Accelerated Failure Time (AFT) model using cooperative learning
#' with two sets of predictors. The function fits a regularized survival model
#' across multiple correlation (rho) values and returns coefficients, predictions,
#' and performance metrics.
#'
#' @param x_training1 A matrix or data frame of the first set of training predictors (e.g., first omic features).
#' @param x_training2 A matrix or data frame of the second set of training predictors (e.g., second omic features).
#' @param y_training A data frame containing survival outcomes with columns:
#'   \itemize{
#'     \item \code{overall_survival}: Numeric survival time
#'     \item \code{status}: Numeric event indicator (1 = event, 0 = censored)
#'   }
#' @param model Character string specifying the AFT distribution. Default is "weibull".
#'   Other options include "lognormal" and "loglogistic".
#' @param rho_values Numeric vector of correlation penalty parameters to evaluate.
#'   Default is \code{c(0.25, 0.5, 0.75)}. Controls the cooperation between the two
#'   predictor sets.
#' @param ncore_max_rho Integer specifying the maximum number of cores to use for
#'   parallel processing across rho values. Default is 10.
#' @param ncore_max_cv Integer specifying the maximum number of cores to use for
#'   cross-validation. Default is 10.
#' @param seed Integer seed for reproducibility. Default is 123.
#'
#' @return A list containing:
#'   \item{final_model_coefs}{Data frame of coefficient estimates for each rho value.
#'     Rows correspond to features from both predictor sets, columns to rho values.}
#'   \item{sigma_est}{Numeric estimated scale parameter (sigma) from the AFT model.}
#'   \item{rho_values}{Numeric vector of rho values used in training (same as input).}
#'   \item{pred}{Training set predictions for each rho value.}
#'   \item{cindex}{Concordance indices (C-index) on training data for each rho value.}
#'   \item{selected_features}{Named list of selected (non-zero) features for each rho value.}
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Prepares survival data (time, event indicator, log-transformed time)
#'   \item Estimates the scale parameter (sigma) using \code{survival::survreg}
#'   \item Fits cooperative AFT models across specified rho values using \code{aft_coop}
#'   \item Generates predictions and calculates C-index for each rho value on training data
#'   \item Identifies selected features (non-zero coefficients) for each rho value
#' }
#'
#' @export

train_AFTCoop = function(x_training1, x_training2, y_training,
                         model = "weibull", rho_values = c(0.25, 0.5, 0.75),
                         ncore_max_rho = 10, ncore_max_cv = 10, seed = 123) {
  set.seed(seed)

  # Surv
  time = as.numeric(y_training$overall_survival)
  event = as.numeric(y_training$status)
  log_time = log(time)
  y_surv = survival::Surv(time, event)

  # check
  x1 = as.matrix(x_training1)
  x2 = as.matrix(x_training2)

  # sigma
  fit_survreg = survival::survreg(Surv(time, event) ~ 1, dist = model, scale = 0)
  sigma_est = exp(fit_survreg$icoef[2])
  message(paste("Sigma:", sigma_est))

  # training
  message("AFTCoop training...")
  final_model_coefs = aft_coop(
    U = x1,
    Z = x2,
    Y = log_time,
    delta = event,
    model = model,
    case = "coop",
    lam_min = TRUE,
    sigma = sigma_est,
    rho_values = rho_values,
    ncore_max_rho = ncore_max_rho,
    ncore_max_cv = ncore_max_cv,
    iplot = FALSE
  )

  rownames(final_model_coefs) = colnames(cbind(x1, x2))
  colnames(final_model_coefs) = paste0("beta_rho_", rho_values)
  final_model_coefs = as.data.frame(final_model_coefs)

  # pred
  pred_train_list = list()
  c_index_train_list = list()
  selected_features_list = list()

  for (i in 1:length(rho_values)) {
    rho_name = colnames(final_model_coefs)[i]
    beta_coefs = final_model_coefs[, i]

    # Predictions
    pred = predict_aft_coop(U = x1, Z = x2, beta = beta_coefs, case = "coop")
    pred_train_list[[rho_name]] = pred

    # C-index
    c_index = survival::concordance(y_surv ~ pred)
    c_index_train_list[[rho_name]] = c_index$concordance

    # Feature selected
    selected_features_list[[rho_name]] = rownames(final_model_coefs)[beta_coefs != 0]
  }

  return(list(
    final_model_coefs = final_model_coefs,
    sigma_est = sigma_est,
    rho_values = rho_values,
    pred = pred_train_list,
    cindex = c_index_train_list,
    selected_features = selected_features_list
  ))
}

#' Test AFT Cooperative Learning Model
#'
#' Evaluates a trained AFT cooperative learning model on test data by generating
#' predictions and calculating concordance indices for each rho value.
#'
#' @param x_test1 A matrix or data frame of the first set of test predictors.
#'   Must have the same features (columns) as \code{x_training1} used in training.
#' @param x_test2 A matrix or data frame of the second set of test predictors.
#'   Must have the same features (columns) as \code{x_training2} used in training.
#' @param y_test A data frame containing test survival outcomes with columns:
#'   \itemize{
#'     \item \code{overall_survival}: Numeric survival time
#'     \item \code{status}: Numeric event indicator (1 = event, 0 = censored)
#'   }
#' @param training_result A list object returned by \code{\link{train_AFTCoop}}
#'   containing the trained model coefficients and parameters.
#'
#' @return A list containing:
#'   \item{pred}{Named list of test set predictions for each rho value.
#'     Higher predicted values indicate higher risk (shorter predicted survival).}
#'   \item{cindex}{Named list of concordance indices (C-index) on test data for each rho value.
#'     Values range from 0 to 1, where 0.5 indicates random predictions and 1.0 indicates
#'     perfect concordance.}
#'
#' @details
#' This function applies the coefficients learned during training to new test data.
#' For each rho value in the trained model:
#' \enumerate{
#'   \item Extracts the corresponding coefficient vector
#'   \item Generates risk predictions using \code{predict_aft_coop}
#'   \item Calculates the concordance index (C-index) comparing predictions to actual outcomes
#' }
#'
#' @export
test_AFTCoop = function(x_test1, x_test2, y_test, training_result) {

  # check
  x1_test = as.matrix(x_test1)
  x2_test = as.matrix(x_test2)

  # surv
  time = as.numeric(y_test$overall_survival)
  event = as.numeric(y_test$status)

  y_test_surv = survival::Surv(time, event)

  # rho extraction
  final_model_coefs = training_result$final_model_coefs
  rho_values = training_result$rho_values

  # pred
  pred_test_list = list()
  c_index_test_list = list()

  for (i in 1:length(rho_values)) {
    rho_name = colnames(final_model_coefs)[i]
    beta_coefs = final_model_coefs[, i]

    # Prediction
    pred = predict_aft_coop(U = x1_test, Z = x2_test, beta = beta_coefs, case = "coop")
    pred_test_list[[rho_name]] = pred

    # C-index
    c_index = survival::concordance(y_test_surv ~ pred)
    c_index_test_list[[rho_name]] = c_index$concordance
  }

  return(list(
    pred = pred_test_list,
    cindex = c_index_test_list
  ))
}
