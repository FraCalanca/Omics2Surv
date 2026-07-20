#' Train AFT Cooperative Learning Model
#'
#' Trains an Accelerated Failure Time cooperative learning model on exactly
#' two omics layers. Called internally by \code{\link{train_j_models}}.
#'
#' @param x_training1 Numeric matrix (samples x features) for the first omics
#'   layer.
#' @param x_training2 Numeric matrix (samples x features) for the second omics
#'   layer.
#' @param y_training Data frame with \code{overall_survival} and \code{status}.
#' @param model Character string specifying the AFT distribution:
#'   \code{"weibull"} (default), \code{"lognormal"}, or \code{"loglogistic"}.
#' @param rho_values Numeric vector of cooperation penalty parameters to
#'   evaluate. Default is \code{c(0.25, 0.5, 0.75)}.
#' @param ncore_max_rho Integer maximum cores for parallelisation across rho
#'   values. Default is 10.
#' @param ncore_max_cv Integer maximum cores for cross-validation. Default
#'   is 10.
#' @param seed Integer seed. Default is 123.
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{final_model_coefs}: Data frame of coefficients
#'       (features x rho values).
#'     \item \code{sigma_est}: Estimated AFT scale parameter.
#'     \item \code{rho_values}: Rho values used.
#'     \item \code{pred}: Named list of training predictions, one per rho.
#'     \item \code{selected_features}: Named list of selected features per rho.
#'     \item \code{mU}, \code{mZ}: Column means of the two training matrices,
#'       used to centre the test set consistently at prediction time.
#'   }
#'
#' @seealso \code{\link{train_j_models}}, \code{\link{test_AFTCoop}}
#'
#' @export
train_AFTCoop <- function(
    x_training1,
    x_training2,
    y_training,
    model         = "weibull",
    rho_values    = c(0.25, 0.5, 0.75),
    ncore_max_rho = 10,
    ncore_max_cv  = 10,
    seed          = 123
) {

  set.seed(seed)

  time     <- as.numeric(y_training$overall_survival)
  event    <- as.numeric(y_training$status)
  log_time <- log(time)

  x1 <- as.matrix(x_training1)
  x2 <- as.matrix(x_training2)

  # ── Estimate sigma ────────────────────────────────────────────────
  fit_survreg <- survival::survreg(
    survival::Surv(time, event) ~ 1,
    dist  = model,
    scale = 0
  )
  sigma_est <- exp(fit_survreg$icoef[[2]])
  message("Sigma: ", round(sigma_est, 4))

  # ── Train ─────────────────────────────────────────────────────────
  message("AFTCoop training...")
  final_model_coefs <- aft_coop(
    U             = x1,
    Z             = x2,
    Y             = log_time,
    delta         = event,
    model         = model,
    case          = "coop",
    lam_min       = TRUE,
    sigma         = sigma_est,
    rho_values    = rho_values,
    ncore_max_rho = ncore_max_rho,
    ncore_max_cv  = ncore_max_cv,
    iplot         = FALSE
  )

  rownames(final_model_coefs) <- colnames(cbind(x1, x2))
  colnames(final_model_coefs) <- paste0("beta_rho_", rho_values)
  final_model_coefs           <- as.data.frame(final_model_coefs)

  # ── Predictions and selected features per rho ─────────────────────
  pred_list     <- vector("list", length(rho_values))
  selected_list <- vector("list", length(rho_values))
  rho_names     <- colnames(final_model_coefs)

  for (i in seq_along(rho_values)) {
    beta_i             <- final_model_coefs[, i]
    pred_list[[i]]     <- predict_aft_coop(U = x1, Z = x2,
                                           beta = beta_i, case = "coop")
    selected_list[[i]] <- rownames(final_model_coefs)[beta_i != 0]
  }

  names(pred_list) <- names(selected_list) <- rho_names

  list(
    final_model_coefs = final_model_coefs,
    sigma_est         = sigma_est,
    rho_values        = rho_values,
    pred              = pred_list,
    selected_features = selected_list,
    mU                = colMeans(x1),   # training means for test centering
    mZ                = colMeans(x2)
  )
}


#' Test AFT Cooperative Learning Model
#'
#' Applies a trained AFTCoop model to test data, generating predictions for
#' each rho value. The test matrices are centred using the training column
#' means (\code{mU}, \code{mZ}) stored by \code{\link{train_AFTCoop}}, so no
#' test-set information leaks into the linear predictor.
#'
#' @param x_test1 Numeric matrix (samples x features) for the first omics
#'   layer. Must match features used in training.
#' @param x_test2 Numeric matrix (samples x features) for the second omics
#'   layer. Must match features used in training.
#' @param training_result List returned by \code{\link{train_AFTCoop}}.
#'
#' @return A list with \code{pred}: named list of numeric prediction vectors,
#'   one per rho value.
#'
#' @seealso \code{\link{train_AFTCoop}}, \code{\link{predict_j_models}}
#'
#' @export
test_AFTCoop <- function(x_test1, x_test2, training_result) {

  x1 <- as.matrix(x_test1)
  x2 <- as.matrix(x_test2)

  final_model_coefs <- training_result$final_model_coefs
  rho_names         <- colnames(final_model_coefs)

  pred_list <- vector("list", length(rho_names))

  for (i in seq_along(rho_names)) {
    beta_i         <- final_model_coefs[, i]
    pred_list[[i]] <- predict_aft_coop(
      U    = x1, Z = x2,
      mU   = training_result$mU,   # centre on TRAINING means
      mZ   = training_result$mZ,
      beta = beta_i,
      case = "coop"
    )
  }

  names(pred_list) <- rho_names

  list(pred = pred_list)
}
