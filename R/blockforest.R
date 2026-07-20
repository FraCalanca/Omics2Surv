#' Train Block Random Forest for Multi-Omics Survival Analysis
#'
#' Trains a block random forest survival model. Feature blocks are derived
#' automatically from the names of \code{X_list}, prefixing each omic's
#' column names. Can be called directly or via \code{\link{train_j_models}}.
#'
#' @param X_list Named list of omics matrices (samples x features). Names are
#'   used to derive block prefixes automatically (e.g. a layer named
#'   \code{"mRNA"} maps to columns starting with \code{"mRNA_"}).
#'   Column names in each matrix must already carry the prefix
#'   \code{"<omic_name>_"} (set automatically by \code{\link{train_j_models}};
#'   must be set manually when calling directly).
#' @param y_training Data frame with \code{overall_survival} and \code{status}.
#' @param threads Integer number of parallel threads. Default is 10.
#' @param always.split.variables Optional list of variable names considered at
#'   every split. Default is \code{NULL}.
#' @param n_top_features Integer number of top important features to return.
#'   Default is 100.
#' @param seed Integer seed. Default is 123.
#' @param ... Additional arguments passed to \code{blockForest::blockForest}.
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{model}: Trained blockForest model.
#'     \item \code{pred}: Numeric vector of risk scores for training samples.
#'     \item \code{selected_features}: Top \code{n_top_features} by impurity.
#'     \item \code{extra}: List with \code{variable_importance} and
#'       \code{blocks}.
#'   }
#'
#' @details
#' Column names are prefixed internally before concatenation so that block
#' membership can be recovered from column names. If calling directly (without
#' \code{\link{train_j_models}}), ensure column names in each matrix of
#' \code{X_list} do not already carry the prefix — they will be added here.
#'
#' @seealso \code{\link{test_blockForest}}, \code{\link{train_j_models}}
#'
#' @export
train_blockForest <- function(
    X_list,
    y_training,
    threads               = 10,
    always.split.variables = NULL,
    n_top_features        = 100,
    seed                  = 123,
    ...
) {

  if (!requireNamespace("blockForest", quietly = TRUE)) {
    stop("Package 'blockForest' is required.")
  }

  set.seed(seed)

  stopifnot(is.list(X_list), !is.null(names(X_list)))

  # ── Prefix column names by omic name ─────────────────────────────
  X_list_prefixed <- lapply(names(X_list), function(nm) {
    x <- X_list[[nm]]
    colnames(x) <- paste0(nm, "_", colnames(x))
    x
  })
  names(X_list_prefixed) <- names(X_list)

  x_training    <- do.call(cbind, X_list_prefixed)
  feature_names <- colnames(x_training)
  prefixes      <- paste0(names(X_list), "_")

  # ── Build blocks ──────────────────────────────────────────────────
  blocks <- lapply(setNames(prefixes, names(X_list)), function(p) {
    idx <- which(startsWith(feature_names, p))
    if (length(idx) == 0L) warning("No features found for prefix: '", p, "'")
    idx
  })

  if (all(sapply(blocks, length) == 0L)) {
    stop("No blocks created. Check names(X_list).")
  }

  time  <- as.numeric(y_training$overall_survival)
  event <- as.numeric(y_training$status)

  if (length(time) != nrow(x_training)) {
    stop("Mismatch between samples in X_list and y_training.")
  }

  if (!is.null(always.split.variables) && !is.list(always.split.variables)) {
    warning("'always.split.variables' must be a list; setting to NULL.")
    always.split.variables <- NULL
  }

  mtry     <- sapply(blocks, function(idx) max(1L, floor(sqrt(length(idx)))))
  BF_train <- data.frame(time = time, status = event, as.data.frame(x_training))
  colnames(BF_train) <- make.names(colnames(BF_train), unique = TRUE)

  # ── Train ─────────────────────────────────────────────────────────
  model_BF <- blockForest::blockForest(
    formula                = survival::Surv(time, status) ~ .,
    data                   = BF_train,
    blocks                 = blocks,
    block.weights          = rep(1, length(blocks)),
    mtry                   = mtry,
    num.threads            = threads,
    importance             = "impurity",
    splitrule              = "logrank",
    write.forest           = TRUE,
    always.split.variables = always.split.variables,
    seed                   = seed,
    ...
  )

  pred           <- model_BF$chf[, ncol(model_BF$chf)]
  names(pred) <- rownames(BF_train)
  var_imp        <- model_BF$variable.importance
  top_features   <- names(head(sort(var_imp, decreasing = TRUE), n_top_features))

  list(
    model             = model_BF,
    pred              = pred,
    selected_features = top_features,
    extra             = list(
      variable_importance = var_imp,
      blocks              = blocks
    )
  )
}


#' Test Block Random Forest Model
#'
#' Applies a trained block random forest model to new data. Can be called
#' directly or via \code{\link{predict_j_models}}.
#'
#' @param X_list Named list of omics matrices (samples x features) for the
#'   test set. Names must match those used in training. Column names will be
#'   prefixed automatically as in training.
#' @param training_result List returned by \code{\link{train_blockForest}}.
#' @param threads Integer number of parallel threads. Default is 10.
#'
#' @return A list with \code{pred}: numeric vector of risk scores.
#'
#' @seealso \code{\link{train_blockForest}}, \code{\link{predict_j_models}}
#'
#' @export
test_blockForest <- function(X_list, training_result, threads = 10) {

  # ── Prefix column names consistently with training ────────────────
  X_list_prefixed <- lapply(names(X_list), function(nm) {
    x <- X_list[[nm]]
    colnames(x) <- paste0(nm, "_", colnames(x))
    x
  })
  names(X_list_prefixed) <- names(X_list)

  x_test  <- do.call(cbind, X_list_prefixed)
  BF_test <- as.data.frame(x_test)
  colnames(BF_test) <- make.names(colnames(BF_test), unique = TRUE)

  pred_BF <- predict(training_result$model,
                     data        = BF_test,
                     num.threads = threads)

  if (is.null(pred_BF$chf)) {
    stop("Prediction object does not contain CHF values.")
  }

  pred <- pred_BF$chf[, ncol(pred_BF$chf)]
  names(pred) <- rownames(BF_test)
  list(pred = pred)
}
