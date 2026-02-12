#' Train Block Random Forest for Survival Analysis
#'
#' Trains a block random forest model for survival data using the blockForest algorithm.
#' Features are organized into blocks based on specified prefixes (e.g., clinical, genomic, transcriptomic),
#' allowing for structured variable selection and importance ranking.
#'
#' @param x_training A matrix or data frame of training predictors.
#' @param y_training A data frame containing survival outcomes with columns:
#'   \itemize{
#'     \item \code{overall_survival}: Numeric survival time
#'     \item \code{status}: Numeric event indicator (1 = event, 0 = censored)
#'   }
#' @param prefixes Character vector of prefixes used to define feature blocks. Each prefix
#'   identifies a group of related variables (e.g., \code{c("clinical_", "omic1_", "omic2_")}).
#'   This parameter is required and must not be NULL.
#' @param threads Integer specifying the number of parallel threads to use. Default is 10.
#' @param always.split.variables Optional list of variable names that should be considered
#'   at every split. Useful for forcing important clinical variables to be evaluated.
#'   Must be a list; set to NULL to disable. Default is NULL.
#' @param n_top_genes Integer specifying the number of top important features to select
#'   based on variable importance. Default is 100.
#' @param seed Integer seed for reproducibility. Default is 8240.
#'
#' @return A list containing:
#'   \item{model}{The trained blockForest model object.}
#'   \item{cindex}{Numeric training concordance index (C-index), calculated as 1 - prediction error.}
#'   \item{pred}{Numeric vector of linear predictors (prognostic scores) for training samples,
#'     extracted from the cumulative hazard function.}
#'   \item{selected_features}{Character vector of the top N most important features based on
#'     impurity-based variable importance, where N = \code{n_top_genes}.}
#'   \item{extra}{List containing additional information:
#'     \itemize{
#'       \item \code{variable_importance}: Named numeric vector of importance scores for all features
#'       \item \code{blocks}: List of feature block indices as defined by prefixes
#'     }}
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Organizes features into blocks based on the specified prefixes
#'   \item Calculates block-specific mtry values using the square root rule: \code{sqrt(block_size)}
#'   \item Trains a survival random forest using the logrank splitting rule
#'   \item Computes variable importance using impurity reduction
#'   \item Extracts training predictions and performance metrics
#'   \item Selects top important features for downstream analysis
#' }
#'
#' @export

train_blockForest = function(x_training, y_training,
                             prefixes = NULL,
                             threads = 10, always.split.variables = NULL, n_top_genes=100,
                             seed = 8240) {

  # Check required packages
  if (!requireNamespace("blockForest", quietly = TRUE)) stop("Package 'blockForest' is required.")
  if (!requireNamespace("survival", quietly = TRUE)) stop("Package 'survival' is required.")

  set.seed(seed)

  # Check the prefix parameter
  if (is.null(prefixes) || !is.character(prefixes)) {
    stop("The prefixes parameter must contain character prefixes, such as variables in data.")
  }

  # Define blocks based on prefixes
  feature_names = colnames(x_training)

  # Create a list where each element is a vector of column indices matching a prefix.
  blocks = sapply(prefixes, function(p) {
    idx = which(startsWith(feature_names, p))
    if (length(idx) == 0) {
      warning(paste("No features found for prefix:", p))
    }
    return(idx)
  }, simplify = FALSE, USE.NAMES = TRUE) # Keep list structure and names

  # Check if blocks were created
  if (length(blocks) == 0 || all(sapply(blocks, length) == 0)) {
    stop("No blocks created. Check prefixes and column names in x_training.")
  }

  # Prepare survival outcome (time and status)
  time = as.numeric(y_training$overall_survival)
  event = as.numeric(y_training$status)
  if (length(time) != nrow(x_training)) stop("Mismatch between number of samples in x_training and y_training.")

  # Combine survival data and predictors into a single data frame
  BF_train = data.frame(time = time, status = event, as.data.frame(x_training))
  colnames(BF_train) = make.names(colnames(BF_train), unique = TRUE)

  # Check always.split.variables parameter
  if (!is.null(always.split.variables) && !is.list(always.split.variables)) {
    message("always.split.variables needs to be a list of the clinical variables present")
    always.split.variables <- NULL
  }

  # Compute mtry for each block (sqrt rule)
  mtry = sapply(blocks, function(idx) max(1, floor(sqrt(length(idx)))))

  # Train the blockForest model
  model_BF = blockForest::blockForest(
    formula = survival::Surv(time, status) ~ .,
    data = BF_train,
    blocks = blocks,
    block.weights = rep(1, length(blocks)),
    mtry = mtry,
    num.threads = threads,
    importance = "impurity",
    splitrule = "logrank",
    write.forest = TRUE,
    always.split.variables = always.split.variables,
    seed = seed
  )

  # Compute training C-index
  if (is.null(model_BF$prediction.error)) {
    stop("Model does not contain 'prediction.error'; cannot compute training C-index.")
  }
  pred_err = model_BF$prediction.error
  train_cindex = 1 - as.numeric(pred_err)

  #Linear predictors (prognostic scores)
  lin_pred_train_BF = model_BF$chf[, ncol(model_BF$chf)]

  #Variable importance and selection of most important genes
  var_imp = model_BF$variable.importance
  var_imp_sorted = sort(var_imp, decreasing = TRUE)
  top_genes = names(head(var_imp_sorted, n_top_genes))

  # Return results
  return(list(
    model = model_BF,
    cindex = train_cindex,
    pred = lin_pred_train_BF,
    selected_features = top_genes,
    extra = list(
      variable_importance = var_imp,
      blocks = blocks
    )
  ))
}

#' Test Block Random Forest Model
#'
#' Evaluates a trained block random forest model on test data by generating risk
#' predictions and calculating the concordance index.
#'
#' @param x_test A matrix or data frame of test predictors. Must have the same
#'   features (columns) as \code{x_training} used in training, including matching
#'   column names and prefixes.
#' @param y_test A data frame containing test survival outcomes with columns:
#'   \itemize{
#'     \item \code{overall_survival}: Numeric survival time
#'     \item \code{status}: Numeric event indicator (1 = event, 0 = censored)
#'   }
#' @param training_result A list object returned by \code{\link{train_blockForest}}
#'   containing the trained model and associated parameters. Note: The function
#'   expects \code{training_result$model_BF}, but \code{train_blockForest} returns
#'   \code{training_result$model}. Ensure compatibility or adjust accordingly.
#' @param num_threads Integer specifying the number of parallel threads to use for
#'   prediction. Default is 10.
#'
#' @return A list containing:
#'   \item{pred}{Numeric vector of linear predictors (risk scores) for test samples,
#'     extracted from the final column of the cumulative hazard function (CHF).
#'     Higher values indicate higher risk (shorter predicted survival).}
#'   \item{cindex}{Numeric concordance index (C-index) on test data. Values range from
#'     0 to 1, where 0.5 indicates random predictions and 1.0 indicates perfect concordance.
#'     Calculated with \code{reverse = TRUE} to account for risk score directionality.}
#'
#' @details
#' This function applies the trained block random forest model to new test data:
#' \enumerate{
#'   \item Prepares test data in the same format as training data
#'   \item Generates predictions using the trained model
#'   \item Extracts risk scores from the cumulative hazard function (CHF)
#'   \item Calculates the concordance index to evaluate prediction accuracy
#' }
#'
#' @export
test_blockForest = function(x_test, y_test, training_result, num_threads = 10) {

  # Prepare survival outcome (time and status)
  time = as.numeric(y_test$overall_survival)
  event = as.numeric(y_test$status)
  if (length(time) != nrow(x_test)) stop("Mismatch between number of samples in x_test and y_test.")

  # Combine survival data and predictors
  BF_test = data.frame(time = time, status = event, as.data.frame(x_test))
  colnames(BF_test) = make.names(colnames(BF_test), unique = TRUE)

  # Predict on the test set
  pred_BF = predict(training_result$model,
                                             data = BF_test,
                                             num.threads = num_threads)

  # Extract risk scores from cumulative hazard function (CHF)
  if (is.null(pred_BF$chf)) stop("Prediction object does not contain CHF values.")
  lp_test = pred_BF$chf[, ncol(pred_BF$chf)]

  # Compute test C-index
  surv_obj_test = survival::Surv(time = BF_test$time, event = BF_test$status)

  c_index_test_obj = survival::concordance(surv_obj_test ~ lp_test, reverse = TRUE)
  c_index_test = c_index_test_obj$concordance

  # Return results
  return(list(
    pred = lp_test,
    cindex = c_index_test
  ))

}
