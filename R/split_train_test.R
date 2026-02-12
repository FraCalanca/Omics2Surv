#' Split Data into Training and Testing Sets
#'
#' Randomly partitions a design matrix and corresponding survival data into
#' training and testing subsets based on a specified proportion.
#'
#' @param design_matrix A matrix or data frame of predictor variables (omics data),
#'   where rows represent samples.
#' @param surv_data A matrix or data frame containing survival outcomes,
#'   with row names matching `design_matrix`.
#' @param seed An integer used to set the random seed for reproducibility.
#'   Default is 4820.
#' @param train_fraction A numeric value between 0 and 1 indicating the
#'   proportion of samples to include in the training set.
#'
#' @return A named list containing four elements:
#' \itemize{
#'   \item \code{train_omic}: Training subset of the design matrix.
#'   \item \code{test_omic}: Testing subset of the design matrix.
#'   \item \code{train_survival}: Training subset of the survival data.
#'   \item \code{test_survival}: Testing subset of the survival data.
#' }
#' @export
split_train_test = function(design_matrix, surv_data, seed = 4820, train_fraction = NULL) {

  set.seed(seed)

  sample_order = rownames(design_matrix)

  surv_data = surv_data[sample_order, , drop = FALSE]

  sample_count = nrow(design_matrix)

  train_samples = sample(seq_len(sample_count), size = round(sample_count * train_fraction), replace = FALSE)
  test_samples = setdiff(seq_len(sample_count), train_samples)

  train_omic = design_matrix[train_samples, ,drop = FALSE]
  test_omic = design_matrix[test_samples, ,drop = FALSE]

  train_survival = surv_data[train_samples, , drop = FALSE]
  test_survival = surv_data[test_samples, , drop = FALSE]

  message("Training set: ", nrow(train_omic), " samples")
  message("Test set: ", nrow(test_omic), " samples")

  return(list(train_omic = train_omic, test_omic = test_omic,
              train_survival = train_survival, test_survival = test_survival))
}
