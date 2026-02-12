#' Train Cooperative Learning Model for Multi-Omics Survival Analysis
#'
#' Trains a cooperative Lasso (CoopLasso) Cox regression model that integrates
#' multiple omics datasets for survival prediction. The method encourages cooperation
#' between different data modalities while performing feature selection through
#' L1 regularization.
#'
#' @param X_list A named list of 2 or 3 omics data matrices or data frames. Each element
#'   represents a different omics layer (e.g., transcriptomics, proteomics, methylation).
#'   All matrices must have:
#'   \itemize{
#'     \item Samples as rows with matching row names across all omics layers
#'     \item Features as columns
#'   }
#' @param y_training A data frame containing survival outcomes with:
#'   \itemize{
#'     \item Row names matching sample identifiers in \code{X_list}
#'     \item Column \code{overall_survival}: Numeric survival time
#'     \item Column \code{status}: Numeric event indicator (1 = event, 0 = censored)
#'   }
#' @param lambdaVector Numeric vector of lambda values to evaluate during cross-validation.
#'   Default is \code{seq(0.1, 0.9, 0.1)}. Lambda controls the strength of regularization.
#' @param folds Integer specifying the number of cross-validation folds for lambda selection.
#'   Default is 5.
#' @param alpha Numeric value controlling the agreement between layers.
#'   \code{alpha = 1} for 2 layers, \code{alpha = 0.5} for 3 layers. Default is 1.
#' @param maxit Integer maximum number of iterations for the optimization algorithm.
#'   Default is 100.
#' @param seed Integer seed for reproducibility of train/validation split and cross-validation.
#'   Default is 123.
#'
#' @return A list containing:
#'   \item{model}{The trained cooplasso model object containing coefficient vector \code{b}
#'     and other model parameters.}
#'   \item{selected_lambda}{Numeric optimal lambda value selected through cross-validation.}
#'   \item{selected_features}{Character vector of feature names with non-zero coefficients
#'     (selected features across all omics layers).}
#'   \item{pred}{Numeric vector of risk scores for training samples. Higher values indicate
#'     higher risk (shorter predicted survival).}
#'   \item{cindex}{Numeric concordance index (C-index) on training data, calculated with
#'     \code{reverse = TRUE} to account for risk score directionality.}
#'   \item{omics}{Character vector of omics layer names from the input \code{X_list}.}
#'
#' @details
#' The function performs the following workflow:
#' \enumerate{
#'   \item Validates input structure (2 or 3 omics layers required)
#'   \item Aligns samples across all omics layers and survival data using common identifiers
#'   \item Splits data into training (90\%) and validation (10\%) sets
#'   \item Performs k-fold cross-validation to select optimal lambda parameter
#'   \item Trains final cooperative Lasso model with selected lambda
#'   \item Computes risk scores and C-index on training data
#'   \item Identifies selected features with non-zero coefficients
#' }
#'
#' @export
train_cooplearning <- function(
    X_list,
    y_training,
    lambdaVector = seq(0.1, 0.9, 0.1),
    folds = 5,
    alpha = 1,
    maxit = 100,
    seed = 123
) {
  set.seed(seed)
  # -----------------------------
  # checks
  # -----------------------------
  stopifnot(is.list(X_list))
  stopifnot(length(X_list) %in% c(2, 3))
  # matching samples
  common_ids <- Reduce(intersect, lapply(X_list, rownames))
  X_list <- lapply(X_list, function(x) x[common_ids, , drop = FALSE])
  y_training <- y_training[common_ids, , drop = FALSE]
  # -----------------------------
  # split train/validation (90/10)
  # -----------------------------
  split <- split_train_test(
    design_matrix = X_list[[1]],
    surv_data = y_training,
    train_fraction = 0.9,
    seed = seed
  )
  train_ids <- rownames(split$train_omic)
  valid_ids <- rownames(split$test_omic)
  X_train_list <- lapply(X_list, function(x) x[train_ids, , drop = FALSE])
  X_valid_list <- lapply(X_list, function(x) x[valid_ids, , drop = FALSE])
  y_train <- split$train_survival
  # -----------------------------
  # survival objects
  # -----------------------------
  time  <- as.numeric(y_train$overall_survival)
  event <- as.numeric(y_train$status)
  # -----------------------------
  # unpack for cooplasso
  # -----------------------------
  Z1 <- X_train_list[[1]]
  Z2 <- X_train_list[[2]]
  Z3 <- if (length(X_train_list) == 3) X_train_list[[3]] else NULL
  Z_valid <- do.call(cbind, X_valid_list)
  # -----------------------------
  # cross-validation lambda
  # -----------------------------
  message("Cross-validation for CoopCox optimal lambda...")
  opt_lambda <- crossval(
    T = time,
    y = event,
    X1 = Z1,
    X2 = Z2,
    X3 = Z3,
    lambdaVector = lambdaVector,
    folds = folds,
    alpha = alpha,
    maxit = maxit
  )
  message("Optimal lambda: ", opt_lambda)
  # -----------------------------
  # training
  # -----------------------------
  message("Training CoopCox...")
  final_model <- cooplasso(
    Z1 = Z1,
    Z2 = Z2,
    Z3 = Z3,
    T = time,
    delta = event,
    validation = Z_valid,
    lambda = opt_lambda,
    alpha = alpha,
    maxit = maxit,
    normalize = TRUE
  )
  # -----------------------------
  # risk score TRAIN
  # -----------------------------
  X_train_combined <- do.call(cbind, X_train_list)
  y_surv_train <- survival::Surv(
    y_training[train_ids, "overall_survival"],
    y_training[train_ids, "status"]
  )
  risk_score_train <- as.numeric(X_train_combined %*% final_model$b)
  train_cindex <- survival::concordance(
    y_surv_train ~ risk_score_train,
    reverse = TRUE
  )$concordance
  # -----------------------------
  # selected features
  # -----------------------------
  coef_vec <- as.numeric(final_model$b)
  names(coef_vec) <- colnames(X_train_combined)
  selected_features <- names(coef_vec)[coef_vec != 0]
  # -----------------------------
  # output
  # -----------------------------
  return(list(
    model = final_model,
    selected_lambda = opt_lambda,
    selected_features = selected_features,
    pred = risk_score_train,
    cindex = train_cindex,
    omics = names(X_list)
  ))
}

#' Test Cooperative Learning Model
#'
#' Evaluates a trained cooperative learning model on test data by generating risk
#' predictions and calculating the concordance index.
#'
#' @param X_list_test A named list of omics data matrices or data frames for testing.
#'   Must contain the same omics layers (in the same order and with matching feature names)
#'   as used in training. Each element should have:
#'   \itemize{
#'     \item Samples as rows with matching row names across all omics layers
#'     \item Features as columns, matching those in the training data
#'   }
#' @param y_test A data frame containing test survival outcomes with:
#'   \itemize{
#'     \item Row names matching sample identifiers in \code{X_list_test}
#'     \item Column \code{overall_survival}: Numeric survival time
#'     \item Column \code{status}: Numeric event indicator (1 = event, 0 = censored)
#'   }
#' @param training_result A list object returned by \code{\link{train_cooplearning}}
#'   containing the trained model and coefficients.
#'
#' @return A list containing:
#'   \item{pred}{Numeric vector of risk scores for test samples. Higher values indicate
#'     higher risk (shorter predicted survival). Calculated as the linear combination
#'     of features weighted by learned coefficients.}
#'   \item{cindex}{Numeric concordance index (C-index) on test data. Values range from
#'     0 to 1, where 0.5 indicates random predictions and 1.0 indicates perfect concordance.
#'     Calculated with \code{reverse = TRUE} to account for risk score directionality.}
#'
#' @details
#' This function applies the trained cooperative learning model to new test data:
#' \enumerate{
#'   \item Validates that \code{X_list_test} is a list
#'   \item Aligns test samples across all omics layers using common identifiers
#'   \item Combines all omics layers into a single feature matrix
#'   \item Ensures feature names match the training model
#'   \item Computes risk scores as a weighted sum: \code{X_test \%*\% coefficients}
#'   \item Calculates C-index to evaluate prediction performance
#' }
#'
#' The C-index measures the model's discriminative ability—the probability that,
#' for a randomly selected pair of patients, the model correctly ranks which patient
#' will experience the event first.
#'
#' @section Data Requirements:
#' \itemize{
#'   \item Test omics layers must match training layers in number, order, and features
#'   \item Only samples present in all test omics layers are used
#'   \item Feature names must be consistent with training data
#' }
#'
#' @seealso \code{\link{train_cooplearning}} for training the model
#'
#' @export
test_cooplearning <- function( X_list_test,
                               y_test,
                               training_result
) {
  stopifnot(is.list(X_list_test))
  common_ids <- Reduce(intersect, lapply(X_list_test, rownames))
  X_list_test <- lapply(X_list_test, function(x) x[common_ids, , drop = FALSE])
  y_test <- y_test[common_ids, , drop = FALSE]
  X_test_combined <- do.call(cbind, X_list_test)
  if (!is.null(names(training_result$model$b))) {
    colnames(X_test_combined) <- names(training_result$model$b)
  }
  # survival
  y_test_surv <- survival::Surv(
    as.numeric(y_test$overall_survival),
    as.numeric(y_test$status)
  )
  # prediction
  pred_test <- as.numeric(X_test_combined %*% training_result$model$b)
  # C-index
  test_cindex <- survival::concordance(
    y_test_surv ~ pred_test,
    reverse = TRUE
  )$concordance
  return(list(
    pred = pred_test,
    cindex = test_cindex
  ))
}
