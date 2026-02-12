configure_flexynesis_python <- function(
    conda_env = "flexynesisenv",
    conda_path = NULL,
    required = TRUE
) {
  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("reticulate is required")

  # 1. If user provides a specific conda binary path, use it
  if (!is.null(conda_path)) {
    options(reticulate.conda_binary = conda_path)
  }

  # 2. Use the environment (reticulate will search standard paths if conda_path is NULL)
  # It can accept either an env name ("flexynesisenv") or a full path
  reticulate::use_condaenv(conda_env, required = required)

  # 3. Validation
  reticulate::py_run_string("
import importlib.util
if importlib.util.find_spec('flexynesis') is None:
    raise ImportError('flexynesis not installed in the active environment')
")
}
