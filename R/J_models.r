#' Train and Test Joint Multi-Omics Models
#'
#' Applies joint (integrative) multi-omics survival analysis approaches.
#' Trains a selected joint model on training data, evaluates it on test data, and
#' stores the results in the MultiAssayExperiment metadata for downstream analysis.
#'
#' @param mae A MultiAssayExperiment object containing the multi-omics data and
#'   clinical information. Results will be stored in \code{metadata(mae)$analysis$cooperative}.
#' @param split A list object returned by a data splitting function, containing:
#'   \itemize{
#'     \item \code{X_list}: List of omics matrices (samples as rows, features as columns)
#'     \item \code{train_survival}: Data frame with training survival outcomes
#'       (columns: \code{overall_survival}, \code{status})
#'     \item \code{test_survival}: Data frame with test survival outcomes
#'   }
#'   Row names of survival data frames are used to identify train/test samples.
#' @param joint_method Character string specifying the joint modeling method to use.
#'   Must be one of:
#'   \itemize{
#'     \item \code{"CoopCox"}: Cooperative Cox regression (supports 2-3 omics layers)
#'     \item \code{"AFTCoop"}: Cooperative Accelerated Failure Time model (requires exactly 2 omics layers)
#'     \item \code{"blockForest"}: Block random forest (supports any number of omics layers)
#'     \item \code{"flexynesis"}: Flexible multi-omics integration via neural networks
#'   }
#' @param omics_list A list of SummarizedExperiment objects or omics data. Required
#'   specifically for the \code{"flexynesis"} method; may be NULL for other methods.
#' @param seed Integer seed for reproducibility. Default is 123. Passed to the
#'   training functions to ensure consistent results.
#' @param ... Additional arguments passed to the specific training function. These
#'   vary by method:
#'   \itemize{
#'     \item \code{CoopCox}: \code{lambdaVector}, \code{folds}, \code{alpha}, \code{maxit}
#'     \item \code{AFTCoop}: \code{model}, \code{rho_values}, \code{ncore_max_rho}, \code{ncore_max_cv}
#'     \item \code{blockForest}: \code{prefixes}, \code{threads}, \code{always.split.variables}, \code{n_top_genes}
#'     \item \code{flexynesis}: method-specific parameters
#'   }
#'
#' @return The input MultiAssayExperiment object with updated metadata. Training and
#'   test results are stored at:
#'   \code{metadata(mae)$analysis$cooperative[[joint_method]]}
#'   which contains:
#'   \itemize{
#'     \item \code{train}: List of training results including model object, predictions,
#'       C-index, and selected features
#'     \item \code{test}: List of test results including predictions and C-index
#'   }
#'   The exact structure of train/test results depends on the specific method used.
#'
#' @details
#' The function performs a unified workflow for joint multi-omics modeling:
#' \enumerate{
#'   \item Extracts training and test sample identifiers from the split object
#'   \item Subsets each omics layer into training and test sets
#'   \item Trains the selected joint model on training data
#'   \item Evaluates the model on held-out test data
#'   \item Stores both training and test results in the MAE metadata
#' }
#'
#' Each joint method has specific requirements:
#' \describe{
#'   \item{CoopCox}{Handles 2-3 omics layers; performs cooperative regularization
#'     with cross-validation for lambda selection}
#'   \item{AFTCoop}{Requires exactly 2 omics layers; models survival time directly
#'     using accelerated failure time framework with cooperative penalties}
#'   \item{blockForest}{Accepts any number of omics layers; organizes features into
#'     blocks for structured random forest modeling}
#'   \item{flexynesis}{Flexible neural network-based integration; requires the
#'     \code{omics_list} parameter}
#' }
#'
#' @section Method-Specific Requirements:
#' \itemize{
#'   \item \strong{AFTCoop}: Will stop with error if \code{length(X_train_list) != 2}
#'   \item \strong{flexynesis}: Requires \code{omics_list} parameter to be specified
#'   \item All methods require consistent sample identifiers across omics layers
#' }
#'
#' @section Result Storage:
#' Results are stored hierarchically in the MAE metadata:
#' \preformatted{
#' metadata(mae)
#'   └─ analysis
#'      └─ cooperative
#'         └─ [joint_method]
#'            ├─ train (training results)
#'            └─ test (test results)
#' }
#'
#' @seealso
#' \code{\link{train_cooplearning}}, \code{\link{test_cooplearning}} for CoopCox
#' \code{\link{train_AFTCoop}}, \code{\link{test_AFTCoop}} for AFTCoop
#' \code{\link{train_blockForest}}, \code{\link{test_blockForest}} for blockForest
#'
#' @export
J_models <- function(
    mae,
    split,
    joint_method = c("CoopCox", "AFTCoop", "blockForest", "flexynesis"),
    omics_list,
    seed = 123,
    ...
) {

  joint_method <- match.arg(joint_method)
  set.seed(seed)

  # 1. Extract split objects
  X_list <- split$X_list
  train_ids <- rownames(split$train_survival)
  test_ids  <- rownames(split$test_survival)

  X_train_list <- lapply(X_list, `[`, train_ids, , drop = FALSE)
  X_test_list  <- lapply(X_list, `[`, test_ids,  , drop = FALSE)

  y_train <- split$train_survival
  y_test  <- split$test_survival

  # 2. Train/Test Logic

  switch(
    joint_method,

    CoopCox = {
      train_res <- train_cooplearning(X_train_list, y_train, seed = seed, ...)
      test_res  <- test_cooplearning(X_test_list, y_test, training_result = train_res)
    },

    AFTCoop = {
      stopifnot(length(X_train_list) == 2)
      train_res <- train_AFTCoop(X_train_list[[1]], X_train_list[[2]], y_train, seed = seed, ...)
      test_res  <- test_AFTCoop(X_test_list[[1]], X_test_list[[2]], y_test, training_result = train_res)
    },

    blockForest = {
      X_train <- do.call(cbind, X_train_list)
      X_test  <- do.call(cbind, X_test_list)
      train_res <- train_blockForest(X_train, y_train, seed = seed, ...)
      test_res  <- test_blockForest(X_test, y_test, training_result = train_res)
    },

    flexynesis = {
      train_res <- train_flexynesis(X_train_list, X_test_list, y_train, y_test, omic_list= omics_list,...)
      test_res  <- test_flexynesis(training_result = train_res, y_test)
    }
  )

  # 3. Store results

  if (is.null(metadata(mae)$analysis$cooperative))
  metadata(mae)$analysis$cooperative <- list()

  metadata(mae)$analysis$cooperative[[joint_method]] <- list(
    train = train_res,
    test = test_res
  )

  return(mae)
}
