# ══════════════════════════════════════════════════════════════════════════════
# Flexynesis
# ──────────────────────────────────────────────────────────────────────────────

# ── Internal: locate a conda/mamba binary in a portable way ─────────────────────
#
#' Locate a conda (or mamba) executable
#'
#' Resolves the path to a conda-family binary without assuming any particular
#' distribution or install location. Resolution order: explicit \code{conda_path},
#' the \code{RETICULATE_CONDA} environment variable, \code{reticulate::conda_binary()},
#' the system \code{PATH}, and finally a probe of common install roots.
#'
#' @param conda_path Optional character path to a conda/mamba binary. If supplied
#'   it takes precedence and must exist.
#'
#' @return A character scalar with the path to the binary.
#' @keywords internal
#' @noRd
.find_conda_binary <- function(conda_path = NULL) {

  if (!is.null(conda_path)) {
    conda_path <- path.expand(conda_path)
    if (file.exists(conda_path)) return(conda_path)
    stop("'conda_path' was supplied but does not exist: ", conda_path)
  }

  env_bin <- Sys.getenv("RETICULATE_CONDA", unset = "")
  if (nzchar(env_bin) && file.exists(env_bin)) return(env_bin)

  rb <- tryCatch(reticulate::conda_binary(), error = function(e) NULL)
  if (!is.null(rb) && nzchar(rb) && file.exists(rb)) return(rb)

  for (exe in c("conda", "mamba", "micromamba")) {
    hit <- Sys.which(exe)
    if (nzchar(hit)) return(unname(hit))
  }

  roots <- c("~/miniforge3", "~/mambaforge", "~/miniconda3", "~/anaconda3",
             "/opt/miniforge3", "/opt/miniconda3", "/opt/anaconda3",
             "/opt/homebrew/Caskroom/miniforge/base",
             "C:/ProgramData/miniconda3", "C:/ProgramData/Anaconda3")
  candidates <- path.expand(c(
    file.path(roots, "bin",      "conda"),
    file.path(roots, "condabin", "conda"),
    file.path(roots, "Scripts",  "conda.exe")
  ))
  hit <- candidates[file.exists(candidates)]
  if (length(hit) > 0L) return(hit[[1]])

  stop("Could not locate a conda/mamba binary. Pass it explicitly via ",
       "'conda_path', or set the RETICULATE_CONDA environment variable.")
}


# ── Internal: build the `conda run` preamble ────────────────────────────────────
#
#' Assemble the leading \code{run} arguments for a conda-family binary
#'
#' Uses \code{--no-capture-output} for conda/mamba so the child process's exit
#' status is propagated and its output streams live (older \code{conda run}
#' buffers everything and can mask a non-zero exit). \code{micromamba run} does
#' not accept that flag, so it is omitted there.
#'
#' @param conda_bin Path to the binary (see \code{.find_conda_binary}).
#' @param conda_env Target environment name.
#' @return Character vector of leading arguments.
#' @keywords internal
#' @noRd
.conda_run_prefix <- function(conda_bin, conda_env) {
  is_micromamba <- grepl("micromamba", basename(conda_bin), ignore.case = TRUE)
  if (is_micromamba) {
    c("run", "-n", conda_env)
  } else {
    c("run", "--no-capture-output", "-n", conda_env)
  }
}


# ── Internal: run a command, streaming output live AND capturing it ─────────────
#
#' Run a command showing output in real time while still capturing it
#'
#' \code{system2(stdout = TRUE)} buffers all output and only returns it when the
#' process ends, so a long flexynesis run looks frozen. This helper instead lets
#' the output stream to the R console as it is produced (so the user sees
#' progress) while also writing it to a log file that is read back afterwards for
#' the sample-removal check. On Unix it pipes through \code{tee}; if that is
#' unavailable it falls back to a plain captured run.
#'
#' @param conda_bin Path to the conda/mamba binary.
#' @param cmd_args Character vector of already-quoted command arguments.
#' @return A list with \code{status} (integer exit code) and \code{log} (the
#'   captured output as a character vector of lines).
#' @keywords internal
#' @noRd
.run_streaming <- function(conda_bin, cmd_args) {
  logfile <- tempfile("flexynesis_log_", fileext = ".log")
  on.exit(unlink(logfile), add = TRUE)

  has_tee <- nzchar(Sys.which("tee")) && .Platform$OS.type == "unix"

  if (has_tee) {
    # Merge stderr into stdout, stream to console, and copy to logfile via tee.
    full_cmd <- paste(
      shQuote(conda_bin), paste(cmd_args, collapse = " "),
      "2>&1 | tee", shQuote(logfile)
    )
    status <- system(full_cmd)               # intern = FALSE -> live console
    log    <- if (file.exists(logfile)) readLines(logfile, warn = FALSE)
    else character(0)
  } else {
    # Fallback: capture only (no live stream), then print.
    log    <- system2(conda_bin, cmd_args, stdout = TRUE, stderr = TRUE)
    status <- attr(log, "status")
    cat(log, sep = "\n")
  }
  if (is.null(status)) status <- 0L
  list(status = as.integer(status), log = log)
}


# ── Internal: verify flexynesis is importable in the target environment ─────────
#
#' Check that flexynesis is installed in a conda environment
#'
#' @param conda_bin Path to the conda/mamba binary.
#' @param conda_env Name of the conda environment containing flexynesis.
#' @return Invisibly \code{TRUE} if the import succeeds; otherwise an error.
#' @keywords internal
#' @noRd
.check_flexynesis_installed <- function(conda_bin, conda_env) {
  # The python snippet must be shQuote()d: system2() does NOT shell-quote its
  # arguments, so an unquoted "import flexynesis" is split on the space and
  # python receives `-c import` (SyntaxError). shQuote keeps it one token.
  res <- system2(
    conda_bin,
    c(.conda_run_prefix(conda_bin, conda_env),
      "python", "-c", shQuote("import flexynesis")),
    stdout = TRUE, stderr = TRUE
  )
  status <- attr(res, "status")
  if (!is.null(status) && status != 0L) {
    stop("flexynesis is not importable in conda env '", conda_env, "':\n",
         paste(res, collapse = "\n"))
  }
  invisible(TRUE)
}


# ── Internal: guard against silent patient loss ─────────────────────────────────
#
#' Assert that sample identifiers are concordant across layers and outcome
#'
#' flexynesis keeps only samples present in every omics layer AND in clin.csv.
#' This check makes any mismatch a hard error *before* export so patients can
#' never be dropped silently by the intersection step.
#'
#' @param X_list Named list of omics matrices (samples x features).
#' @param y Outcome data frame; row names are sample identifiers.
#' @param cohort Label used in messages ("train" / "test").
#' @keywords internal
#' @noRd
.assert_sample_concordance <- function(X_list, y, cohort) {
  ref <- rownames(y)
  if (is.null(ref) || anyNA(ref) || any(!nzchar(ref))) {
    stop("[", cohort, "] y must have non-empty, non-NA row names (sample IDs).")
  }
  if (anyDuplicated(ref)) {
    stop("[", cohort, "] duplicated sample IDs in y: ",
         paste(unique(ref[duplicated(ref)]), collapse = ", "))
  }
  for (nm in names(X_list)) {
    s <- rownames(X_list[[nm]])
    if (is.null(s)) {
      stop("[", cohort, "] omics layer '", nm,
           "' has no row names (sample IDs).")
    }
    missing_in_layer <- setdiff(ref, s)
    missing_in_y     <- setdiff(s, ref)
    if (length(missing_in_layer) > 0L || length(missing_in_y) > 0L) {
      stop("[", cohort, "] sample mismatch between y and layer '", nm, "'.\n",
           if (length(missing_in_layer))
             paste0("  In y but not in '", nm, "': ",
                    paste(utils::head(missing_in_layer, 10), collapse = ", "),
                    if (length(missing_in_layer) > 10) " ...", "\n"),
           if (length(missing_in_y))
             paste0("  In '", nm, "' but not in y: ",
                    paste(utils::head(missing_in_y, 10), collapse = ", "),
                    if (length(missing_in_y) > 10) " ..."),
           "\nflexynesis would intersect these and drop the offending patients. ",
           "Align the cohorts before calling.")
    }
  }
  invisible(TRUE)
}


# ── Internal: detect dropped samples in flexynesis stdout ───────────────────────
#
#' Scan flexynesis log output for sample-removal notices
#'
#' Raises an error if any layer reports a non-zero "Removed N samples" so that
#' patient loss inside flexynesis surfaces immediately.
#'
#' @param log_lines Character vector of captured flexynesis output.
#' @keywords internal
#' @noRd
.assert_no_samples_removed <- function(log_lines) {
  hits <- grep("Removed\\s+([0-9]+)\\s+samples", log_lines, value = TRUE)
  n <- suppressWarnings(as.integer(sub(".*Removed\\s+([0-9]+)\\s+samples.*",
                                       "\\1", hits)))
  if (any(n > 0L, na.rm = TRUE)) {
    stop("flexynesis removed samples during preprocessing:\n",
         paste(hits[n > 0L], collapse = "\n"),
         "\nThis means patients were dropped. Inspect your input matrices ",
         "(constant/all-NA rows, ID mismatches) before re-running.")
  }
  invisible(TRUE)
}


#' Assert that no sample is constant within a layer
#'
#' flexynesis drops any sample whose values have zero (or NaN) standard deviation
#' across features in a layer (the hard-coded \code{common_mask} in
#' \code{cleanup_data}); this is not configurable. Such a sample is "constant" in
#' that layer - e.g. a copy-neutral genome giving an all-equal CNV row. This
#' check reproduces that mask in R so the offending patient IDs surface before
#' export instead of as an opaque "Removed N samples" line.
#'
#' @param X_list Named list of omics matrices (samples x features).
#' @param cohort Label used in messages.
#' @keywords internal
#' @noRd
.assert_informative_samples <- function(X_list, cohort) {
  for (nm in names(X_list)) {
    m  <- X_list[[nm]]
    sd_per_sample <- apply(m, 1L, function(r) stats::sd(r, na.rm = TRUE))
    bad <- rownames(m)[!is.finite(sd_per_sample) | sd_per_sample == 0]
    if (length(bad) > 0L) {
      stop("[", cohort, "] layer '", nm, "': ", length(bad),
           " sample(s) are constant across features (zero/NaN variance). ",
           "flexynesis would silently drop them. Offending IDs: ",
           paste(utils::head(bad, 15), collapse = ", "),
           if (length(bad) > 15) paste0(" ... (+", length(bad) - 15, " more)"),
           "\nThese carry no information in '", nm, "' (e.g. a flat CNV ",
           "profile). Resolve upstream - drop the layer for them, fix the ",
           "matrix, or exclude them deliberately - before re-running.")
    }
  }
  invisible(TRUE)
}
#'
#' Writes one cohort to disk in the directory structure flexynesis expects.
#'
#' Two layouts are supported:
#' \itemize{
#'   \item \code{"train_test"}: a \code{train/} and a \code{test/} sub-folder,
#'     each with one CSV per omics layer plus \code{clin.csv}. This is the layout
#'     for a normal training run (\code{--data_path}).
#'   \item \code{"flat"}: \code{clin.csv} and the omics CSVs written directly in
#'     \code{outdir}, with no sub-folders. This is the layout required by
#'     inference mode (\code{--data_path_test}).
#' }
#'
#' Omics matrices are transposed on write (features as rows, samples as columns),
#' matching the flexynesis input convention. Inputs are therefore expected to be
#' samples-by-features. Sample concordance across layers and \code{y} is asserted
#' first to prevent silent patient loss.
#'
#' @param X_list Named list of omics matrices (samples x features) for this cohort.
#' @param y Data frame of clinical/survival outcome; row names are sample IDs.
#' @param outdir Output directory (created recursively if needed).
#' @param layout One of \code{"train_test"} or \code{"flat"}.
#' @param X_list2,y2 Only for \code{layout = "train_test"}: the second cohort
#'   written into \code{test/} (the first goes into \code{train/}).
#'
#' @return Invisibly the \code{outdir} path.
#' @export
export_flexynesis_data <- function(
    X_list,
    y,
    outdir,
    layout = c("train_test", "flat"),
    X_list2 = NULL,
    y2      = NULL
) {
  layout <- match.arg(layout)

  .write_cohort <- function(Xl, yy, dir, cohort) {
    .assert_sample_concordance(Xl, yy, cohort)
    .assert_informative_samples(Xl, cohort)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    # clin.csv: sample IDs are the (unnamed) first column; pandas reads it as the
    # index via index_col=0.
    utils::write.csv(yy, file.path(dir, "clin.csv"))
    for (omic in names(Xl)) {
      # Reorder layer rows to match clin order, then transpose to features x samples.
      mat <- Xl[[omic]][rownames(yy), , drop = FALSE]
      utils::write.csv(t(mat), file.path(dir, paste0(omic, ".csv")))
    }
  }

  if (layout == "train_test") {
    if (is.null(X_list2) || is.null(y2)) {
      stop("layout = 'train_test' requires X_list2 and y2 for the test folder.")
    }
    .write_cohort(X_list,  y,  file.path(outdir, "train"), "train")
    .write_cohort(X_list2, y2, file.path(outdir, "test"),  "test")
  } else {
    .write_cohort(X_list, y, outdir, "test")
  }

  invisible(outdir)
}


#' Run the flexynesis CLI (training mode)
#'
#' Builds and executes a flexynesis training command via \code{conda run}, with
#' all internal feature selection disabled by default so pre-filtered data passes
#' through untouched. Output is streamed live; on completion the function checks
#' the log for dropped samples and verifies the expected output files exist.
#'
#' @param data_path Folder holding \code{train/} and \code{test/} (see
#'   \code{\link{export_flexynesis_data}} with \code{layout = "train_test"}).
#' @param omic_list Character vector of omics layer names; becomes
#'   \code{--data_types} and must match the CSV file names.
#' @param surv_time,surv_event Column names in \code{clin.csv} for survival time
#'   and event indicator.
#' @param outdir Output folder for flexynesis results.
#' @param model_class flexynesis model class. Default \code{"DirectPred"}. With
#'   \code{surv_time}/\code{surv_event} set, the DL models (DirectPred,
#'   supervised_vae, MultiTripletNetwork, CrossModalPred, GNN) attach an MLP head
#'   trained with Cox-PH loss, so they predict survival risk scores directly.
#' @param prefix Job prefix for output files. Default \code{"job"}.
#' @param threads Number of CPU threads. Default 10.
#' @param device \code{"cpu"}, \code{"mps"}, \code{"cuda"} or \code{"auto"}.
#' @param use_cv Logical; request 5-fold cross-validation.
#' @param variance_threshold,correlation_threshold,features_top_percentile,features_min
#'   Feature-selection controls. Defaults (0, 1, 100, 0) keep every feature, i.e.
#'   no selection. \code{features_top_percentile = 100} (not 0) is used on
#'   purpose: 0 makes older flexynesis skip the selection log and crash with
#'   \code{KeyError: 'select_features'}, whereas 100 keeps all features via the
#'   Laplacian early-return. Override only if you want flexynesis to filter.
#' @param subsample Training downsample size; 0 (default) keeps every patient.
#' @param finetuning_samples Test samples used for fine-tuning; 0 (default) off.
#' @param feature_importance_method Captum explainer for marker discovery:
#'   \code{"IntegratedGradients"} (default), \code{"GradientShap"} or
#'   \code{"Both"}. Marker finding is left enabled so the importance CSV(s) are
#'   written.
#' @param safetensors Logical; save as SafeTensors + JSON instead of .pth + .joblib.
#' @param extra_args Optional character vector of additional raw CLI tokens.
#' @param conda_env,conda_path Conda environment / binary.
#'
#' @return Invisibly a list of output paths: \code{model}, \code{artifacts},
#'   \code{predicted_labels}, \code{feature_importance} (character vector of the
#'   importance CSV(s) actually written).
#' @export
run_flexynesis <- function(
    data_path,
    omic_list,
    surv_time,
    surv_event,
    outdir,
    model_class             = "DirectPred",
    prefix                  = "job",
    threads                 = 10,
    device                  = "cpu",
    use_cv                  = TRUE,
    variance_threshold      = 0,
    correlation_threshold   = 1,
    features_top_percentile = 100,
    features_min            = 0,
    subsample               = 0,
    finetuning_samples      = 0,
    feature_importance_method = c("IntegratedGradients", "GradientShap", "Both"),
    safetensors             = FALSE,
    extra_args              = character(0),
    conda_env               = "flexynesisenv",
    conda_path              = NULL
) {

  feature_importance_method <- match.arg(feature_importance_method)

  conda_bin <- .find_conda_binary(conda_path)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  cmd_args <- c(
    .conda_run_prefix(conda_bin, conda_env), "flexynesis",
    "--data_path",               shQuote(data_path),
    "--model_class",             model_class,
    "--surv_time_var",           surv_time,
    "--surv_event_var",          surv_event,
    "--data_types",              paste(omic_list, collapse = ","),
    "--threads",                 as.character(threads),
    "--device",                  device,
    "--outdir",                  shQuote(outdir),
    "--prefix",                  prefix,
    # ── feature selection fully OFF by default (see file header) ──
    "--variance_threshold",      as.character(variance_threshold),
    "--correlation_threshold",   as.character(correlation_threshold),
    "--features_top_percentile", as.character(features_top_percentile),
    "--features_min",            as.character(features_min),
    "--subsample",               as.character(subsample),
    "--finetuning_samples",      as.character(finetuning_samples),
    # marker finding stays ENABLED (no --disable_marker_finding) so that
    # <prefix>.feature_importance.<explainer>.csv is written for each target.
    "--feature_importance_method", feature_importance_method
  )
  if (isTRUE(use_cv))      cmd_args <- c(cmd_args, "--use_cv")
  if (isTRUE(safetensors)) cmd_args <- c(cmd_args, "--safetensors")
  if (length(extra_args))  cmd_args <- c(cmd_args, extra_args)

  message("Running flexynesis (training):")
  message(paste(c(conda_bin, cmd_args), collapse = " "))

  run    <- .run_streaming(conda_bin, cmd_args)
  res    <- run$log
  status <- run$status
  if (status != 0L) {
    stop("flexynesis exited with status ", status,
         " - see the Python output above.")
  }

  # Patient-loss guard: fail if flexynesis reports any removed samples.
  .assert_no_samples_removed(res)

  ext_model <- if (safetensors) "final_model.safetensors" else "final_model.pth"
  ext_art   <- if (safetensors) "artifacts.json"          else "artifacts.joblib"
  explainers <- if (feature_importance_method == "Both")
    c("IntegratedGradients", "GradientShap") else feature_importance_method
  importance_paths <- file.path(
    outdir, paste0(prefix, ".feature_importance.", explainers, ".csv")
  )
  paths <- list(
    model            = file.path(outdir, paste0(prefix, ".", ext_model)),
    artifacts        = file.path(outdir, paste0(prefix, ".", ext_art)),
    predicted_labels = file.path(outdir, paste0(prefix, ".predicted_labels.csv")),
    feature_importance = importance_paths[file.exists(importance_paths)]
  )
  missing <- vapply(paths[c("model", "artifacts", "predicted_labels")],
                    function(p) !file.exists(p), logical(1))
  if (any(missing)) {
    stop("flexynesis finished but expected output(s) are missing:\n",
         paste(unlist(paths[missing]), collapse = "\n"),
         "\nCheck the Python output above.")
  }

  invisible(paths)
}


#' Run the flexynesis CLI (inference mode)
#'
#' Predicts on a test-only cohort using a model + artifacts produced by
#' \code{\link{run_flexynesis}}. Inference forces CPU inside flexynesis and reuses
#' the exact features and scalers stored in the artifacts, so the feature
#' selection flags do not apply here.
#'
#' @param pretrained_model Path to \code{.final_model.pth} / \code{.safetensors}.
#' @param artifacts Path to \code{.artifacts.joblib} / \code{.artifacts.json}.
#' @param data_path_test FLAT folder with \code{clin.csv} + omics CSVs (see
#'   \code{\link{export_flexynesis_data}} with \code{layout = "flat"}).
#' @param outdir Output folder for predictions.
#' @param prefix Job prefix. Default \code{"job"}.
#' @param surv_time,surv_event Survival column names. These MUST be passed: in
#'   flexynesis the predicted-labels block is gated by
#'   \code{if any([target_variables, surv_event_var])}, so without
#'   \code{--surv_event_var} inference stops after writing embeddings and never
#'   produces \code{predicted_labels.csv}.
#' @param safetensors Logical; set if the model was saved as SafeTensors.
#' @param extra_args Optional additional raw CLI tokens.
#' @param conda_env,conda_path Conda environment / binary.
#'
#' @return Invisibly the path to \code{<prefix>.predicted_labels.csv}.
#' @export
run_flexynesis_inference <- function(
    pretrained_model,
    artifacts,
    data_path_test,
    outdir,
    surv_time   = "overall_survival",
    surv_event  = "status",
    prefix      = "job",
    safetensors = FALSE,
    extra_args  = character(0),
    conda_env   = "flexynesisenv",
    conda_path  = NULL
) {

  conda_bin <- .find_conda_binary(conda_path)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  for (p in c(pretrained_model, artifacts)) {
    if (!file.exists(p)) stop("Required inference input not found: ", p)
  }

  cmd_args <- c(
    .conda_run_prefix(conda_bin, conda_env), "flexynesis",
    "--pretrained_model", shQuote(pretrained_model),
    "--artifacts",        shQuote(artifacts),
    "--data_path_test",   shQuote(data_path_test),
    "--outdir",           shQuote(outdir),
    "--prefix",           prefix,
    # Required: gates flexynesis's prediction step. Without --surv_event_var the
    # run stops after embeddings and never writes predicted_labels.csv.
    "--surv_time_var",    surv_time,
    "--surv_event_var",   surv_event
  )
  if (isTRUE(safetensors)) cmd_args <- c(cmd_args, "--safetensors")
  if (length(extra_args))  cmd_args <- c(cmd_args, extra_args)

  message("Running flexynesis (inference):")
  message(paste(c(conda_bin, cmd_args), collapse = " "))

  run    <- .run_streaming(conda_bin, cmd_args)
  res    <- run$log
  status <- run$status
  if (status != 0L) {
    stop("flexynesis inference exited with status ", status,
         " - see the Python output above.")
  }
  .assert_no_samples_removed(res)

  csv_path <- file.path(outdir, paste0(prefix, ".predicted_labels.csv"))
  if (!file.exists(csv_path)) {
    stop("flexynesis inference finished but the output file was not found:\n",
         csv_path, "\nCheck the Python output above.")
  }
  invisible(csv_path)
}


#' Import flexynesis predictions
#'
#' Reads a \code{predicted_labels.csv} and returns one prediction per sample for
#' a single split. The file holds one row per (sample, variable) for numeric
#' targets and per (sample, variable, class) for categorical targets, so a naive
#' \code{row.names = sample_id} (as in the original code) would hit duplicate row
#' names. This version filters to a single \code{variable}, de-duplicates by
#' sample, and (optionally) verifies that no expected patient is missing.
#'
#' @param csv_path Path to a flexynesis \code{predicted_labels.csv}.
#' @param split Which split to extract: \code{"test"} (default) or \code{"train"}.
#' @param variable Optional outcome/variable name to filter on. If \code{NULL}
#'   and several are present, the first is used (with a warning).
#' @param expected_samples Optional character vector of sample IDs that must all
#'   be present; missing IDs raise an error (patient-loss guard).
#' @param numeric_pred Logical; coerce predictions to numeric (TRUE for survival
#'   risk scores / regression). Default \code{TRUE}.
#'
#' @return A one-column data frame of predictions (\code{pred}) indexed by sample
#'   identifier, ordered to match \code{expected_samples} when supplied.
#' @export
import_flexynesis_results <- function(
    csv_path,
    split            = c("test", "train"),
    variable         = NULL,
    expected_samples = NULL,
    numeric_pred     = TRUE
) {
  split <- match.arg(split)
  preds <- readr::read_csv(csv_path, show_col_types = FALSE)

  preds <- preds[preds$split == split, , drop = FALSE]
  if (nrow(preds) == 0L) {
    stop("No rows with split == '", split, "' in ", csv_path)
  }

  vars <- unique(preds$variable)
  if (is.null(variable)) {
    variable <- vars[[1]]
    if (length(vars) > 1L) {
      warning("Multiple variables in predictions (",
              paste(vars, collapse = ", "), "); using '", variable, "'.")
    }
  } else if (!variable %in% vars) {
    stop("Requested variable '", variable, "' not found. Available: ",
         paste(vars, collapse = ", "))
  }
  preds <- preds[preds$variable == variable, , drop = FALSE]

  # One row per sample (categorical targets repeat the predicted_label per class).
  preds <- preds[!duplicated(preds$sample_id), , drop = FALSE]

  pred_vals <- preds$predicted_label
  if (isTRUE(numeric_pred)) pred_vals <- as.numeric(pred_vals)
  out <- data.frame(pred = pred_vals, row.names = preds$sample_id)

  if (!is.null(expected_samples)) {
    missing <- setdiff(expected_samples, rownames(out))
    if (length(missing) > 0L) {
      stop("flexynesis did not return predictions for ", length(missing),
           " expected patient(s): ",
           paste(utils::head(missing, 10), collapse = ", "),
           if (length(missing) > 10) " ...",
           "\nThey were dropped during preprocessing.")
    }
    out <- out[expected_samples, , drop = FALSE]
  }
  out
}


#' Import flexynesis feature-importance scores
#'
#' Reads the \code{<prefix>.feature_importance.<explainer>.csv} file(s) produced
#' during training (one per explainer) and stacks them. Columns:
#' \code{target_variable}, \code{target_class}, \code{target_class_label},
#' \code{layer}, \code{name} (feature), \code{importance}, \code{explainer}.
#'
#' @param paths Character vector of feature-importance CSV paths (e.g.
#'   \code{run_flexynesis()$feature_importance}).
#' @return A single data frame, or \code{NULL} if no files exist.
#' @export
import_flexynesis_importance <- function(paths) {
  paths <- paths[file.exists(paths)]
  if (length(paths) == 0L) return(NULL)
  do.call(rbind, lapply(paths, function(p) {
    as.data.frame(readr::read_csv(p, show_col_types = FALSE))
  }))
}


#' Train a flexynesis model on the training cohort only
#'
#' Exports the training cohort, runs flexynesis in training mode with all feature
#' selection disabled, and persists the model and artifacts for later inference.
#' The test cohort is never seen here. flexynesis still requires a \code{test/}
#' folder to start, so a copy of the training cohort is written as a harmless
#' placeholder (identical features => the harmonize step is a no-op and nothing
#' from a real test set can leak into training).
#'
#' @param X_train_list Named list of training omics matrices (samples x features).
#' @param y_train Data frame of training survival data; row names are sample IDs.
#' @param omic_list Character vector of omics layer names to use.
#' @param surv_time,surv_event Survival column names. Defaults
#'   \code{"overall_survival"} / \code{"status"}.
#' @param tmpdir Working directory; a fresh per-call subdirectory by default.
#' @param conda_env,conda_path Conda environment / binary.
#' @param check_install Logical; verify flexynesis importable first (default TRUE).
#' @param safetensors Logical; passed through to \code{\link{run_flexynesis}}.
#' @param ... Additional arguments forwarded to \code{\link{run_flexynesis}}
#'   (e.g. \code{model_class}, \code{device}, \code{use_cv}, \code{threads},
#'   or the feature-selection overrides).
#'
#' @return A list with \code{model} (paths to the saved model/artifacts and the
#'   working dirs), \code{pred} (training-split predictions),
#'   \code{feature_importance} (stacked importance table, or \code{NULL}),
#'   \code{selected_features} (alias of \code{feature_importance}), and
#'   \code{params}.
#'
#' @seealso \code{\link{test_flexynesis}}, \code{\link{run_flexynesis}}
#' @export
train_flexynesis <- function(
    X_train_list,
    y_train,
    omic_list,
    surv_time     = "overall_survival",
    surv_event    = "status",
    tmpdir        = tempfile("flexynesis_"),
    conda_env     = "flexynesisenv",
    conda_path    = NULL,
    check_install = TRUE,
    safetensors   = FALSE,
    ...
) {

  if (isTRUE(check_install)) {
    .check_flexynesis_installed(.find_conda_binary(conda_path), conda_env)
  }

  dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)

  # Export train -> train/, and a copy of train -> test/ (placeholder only).
  export_flexynesis_data(
    X_list  = X_train_list, y  = y_train,
    X_list2 = X_train_list, y2 = y_train,
    outdir  = tmpdir, layout = "train_test"
  )

  results_dir <- file.path(tmpdir, "results")
  out_paths <- run_flexynesis(
    data_path  = tmpdir,
    omic_list  = omic_list,
    surv_time  = surv_time,
    surv_event = surv_event,
    outdir     = results_dir,
    conda_env  = conda_env,
    conda_path = conda_path,
    safetensors = safetensors,
    ...
  )

  pred_train <- import_flexynesis_results(
    out_paths$predicted_labels,
    split            = "train",
    variable         = surv_event,
    expected_samples = rownames(y_train),
    numeric_pred     = TRUE
  )

  feature_importance <- import_flexynesis_importance(out_paths$feature_importance)

  list(
    model = list(
      pretrained_model = out_paths$model,
      artifacts        = out_paths$artifacts,
      results_dir      = results_dir,
      tmpdir           = tmpdir,
      safetensors      = safetensors,
      omic_list        = omic_list,
      surv_event       = surv_event,
      surv_time        = surv_time,
      conda_env        = conda_env,
      conda_path       = conda_path
    ),
    pred               = pred_train,
    feature_importance = feature_importance,
    selected_features  = feature_importance,
    params             = list(...)
  )
}


#' Predict with a trained flexynesis model on the test cohort
#'
#' Exports the test cohort to the flat inference layout and runs flexynesis in
#' inference mode against the model + artifacts from \code{\link{train_flexynesis}}.
#'
#' @param training_result List returned by \code{\link{train_flexynesis}}.
#' @param X_test_list Named list of test omics matrices (samples x features).
#' @param y_test Data frame of test survival data; row names are sample IDs.
#'   \code{clin.csv} is written from this so known labels are available, but it is
#'   not required for prediction itself.
#' @param tmpdir Working directory for the exported test data; fresh by default.
#'
#' @return A list with \code{pred}: the test-set predictions, ordered to match
#'   \code{rownames(y_test)}.
#' @seealso \code{\link{train_flexynesis}}, \code{\link{run_flexynesis_inference}}
#' @export
test_flexynesis <- function(
    training_result,
    X_test_list,
    y_test,
    tmpdir = tempfile("flexynesis_infer_")
) {
  m <- training_result$model
  if (is.null(m$pretrained_model) || is.null(m$artifacts)) {
    stop("training_result$model is missing model/artifacts paths. ",
         "Run train_flexynesis() first.")
  }

  infer_dir <- file.path(tmpdir, "test_data")
  dir.create(infer_dir, recursive = TRUE, showWarnings = FALSE)

  # Flat layout for --data_path_test.
  export_flexynesis_data(
    X_list = X_test_list, y = y_test,
    outdir = infer_dir, layout = "flat"
  )

  out_dir  <- file.path(tmpdir, "results")
  csv_path <- run_flexynesis_inference(
    pretrained_model = m$pretrained_model,
    artifacts        = m$artifacts,
    data_path_test   = infer_dir,
    outdir           = out_dir,
    surv_time        = m$surv_time,
    surv_event       = m$surv_event,
    safetensors      = isTRUE(m$safetensors),
    conda_env        = m$conda_env,
    conda_path       = m$conda_path
  )

  pred_test <- import_flexynesis_results(
    csv_path,
    split            = "test",
    variable         = m$surv_event,
    expected_samples = rownames(y_test),
    numeric_pred     = TRUE
  )

  list(pred = pred_test)
}


# ══════════════════════════════════════════════════════════════════════════════
# MAE-level wrappers
# ──────────────────────────────────────────────────────────────────────────────
# These mirror train_j_models() / predict_j_models(): they take a
# MultiAssayExperiment, pull X/y out of it, run the base learners above, compute
# the SAME metrics as the joint methods (Harrell on train; Harrell + Uno on
# test) by reusing the shared internal helpers, and store everything under
# metadata(mae)$analysis$flexynesis.
#
# Shared helpers reused from the joint-models file (same package namespace):
#   .extract_x_list(mae, omic_list)   -> named list of (samples x features)
#   .extract_xy(mae, omic)$y          -> outcome df (overall_survival, status)
#   .compute_train_metrics(train_res, y, method)            -> harrell_train
#   .compute_test_metrics(train_res, test_res, y_tr, y_te, method)
#                                       -> harrell_test + uno_test
# flexynesis (DirectPred + Cox-PH) yields a risk score, so it follows the
# non-AFT branch of those helpers; we pass method = "flexynesis".
# ══════════════════════════════════════════════════════════════════════════════

# Internal: turn a flexynesis pred data frame (pred column, sample row names)
# into a numeric vector aligned to the row order of `y`.
.flexynesis_pred_vec <- function(pred_df, y) {
  as.numeric(pred_df[rownames(y), "pred"])
}


#' Train a flexynesis model from a training MAE
#'
#' MAE-level entry point: extracts the omics layers and survival outcome from
#' \code{mae_train}, trains flexynesis on the training cohort only (all internal
#' feature selection disabled), computes the resubstitution Harrell C-index, and
#' stores the fitted model, training predictions, feature importance and metrics
#' in \code{metadata(mae_train)$analysis$flexynesis}.
#'
#' @param mae_train A \code{MultiAssayExperiment}; \code{colData} must contain
#'   \code{"overall_survival"} and \code{"status"}.
#' @param omic_list Character vector of assay names to use.
#' @param surv_time,surv_event Survival column names. Defaults
#'   \code{"overall_survival"} / \code{"status"}.
#' @param ... Forwarded to \code{\link{train_flexynesis}} (e.g.
#'   \code{model_class}, \code{device}, \code{use_cv}, \code{threads},
#'   \code{conda_env}, feature-selection overrides).
#'
#' @return \code{mae_train} with results at
#'   \code{metadata(mae_train)$analysis$flexynesis}, containing \code{train}
#'   (the \code{\link{train_flexynesis}} result), \code{metrics}
#'   (\code{harrell_train}) and \code{feature_importance}.
#'
#' @seealso \code{\link{predict_flexynesis_models}}, \code{\link{train_flexynesis}}
#' @export
train_flexynesis_models <- function(
    mae_train,
    omic_list,
    surv_time  = "overall_survival",
    surv_event = "status",
    ...
) {

  missing_assays <- setdiff(omic_list, names(mae_train))
  if (length(missing_assays) > 0L) {
    stop("Assay(s) not found in mae_train: ",
         paste(missing_assays, collapse = ", "))
  }

  X_list <- .extract_x_list(mae_train, omic_list)
  y      <- .extract_xy(mae_train, omic_list[[1]])$y

  train_res <- train_flexynesis(
    X_train_list = X_list,
    y_train      = y,
    omic_list    = omic_list,
    surv_time    = surv_time,
    surv_event   = surv_event,
    ...
  )

  # Harrell (resubstitution), via the shared helper. flexynesis is non-AFT.
  metrics <- .compute_train_metrics(
    list(pred = .flexynesis_pred_vec(train_res$pred, y)),
    y, "flexynesis"
  )

  tag <- paste(sort(omic_list), collapse = "_")

  md <- MultiAssayExperiment::metadata(mae_train)
  if (is.null(md$analysis))            md$analysis <- list()
  if (is.null(md$analysis$flexynesis)) md$analysis$flexynesis <- list()
  md$analysis$flexynesis[[tag]] <- list(
    train              = train_res,
    metrics            = metrics,
    feature_importance = train_res$feature_importance
  )
  MultiAssayExperiment::metadata(mae_train) <- md

  mae_train
}


#' Predict from a trained flexynesis model on a new MAE
#'
#' Retrieves the model fitted by \code{\link{train_flexynesis_models}} from
#' \code{mae_train}, runs flexynesis inference on \code{mae_pred}, and computes
#' the Harrell and Uno C-indices (the latter using the training cohort for IPCW,
#' exactly as the joint methods do). Results are stored in
#' \code{metadata(mae_pred)$analysis$flexynesis}.
#'
#' @param mae_pred A \code{MultiAssayExperiment} to predict on.
#' @param mae_train A \code{MultiAssayExperiment} carrying the fitted model in
#'   \code{metadata(mae_train)$analysis$flexynesis}.
#' @param omic_list Character vector of assay names; must match training.
#' @param ... Forwarded to \code{\link{test_flexynesis}}.
#'
#' @return \code{mae_pred} with results at
#'   \code{metadata(mae_pred)$analysis$flexynesis}, containing \code{test} and
#'   \code{metrics} (\code{harrell_test} + \code{uno_test}).
#'
#' @seealso \code{\link{train_flexynesis_models}}, \code{\link{test_flexynesis}}
#' @export
predict_flexynesis_models <- function(
    mae_pred,
    mae_train,
    omic_list,
    ...
) {

  tag <- paste(sort(omic_list), collapse = "_")

  train_res <- MultiAssayExperiment::metadata(mae_train)$analysis$flexynesis[[tag]]$train
  if (is.null(train_res)) {
    stop("No fitted flexynesis model found in mae_train for omics '", tag,
         "'. Run train_flexynesis_models() with the same omic_list first.")
  }

  X_test <- .extract_x_list(mae_pred, omic_list)
  y_test  <- .extract_xy(mae_pred,  omic_list[[1]])$y
  y_train <- .extract_xy(mae_train, omic_list[[1]])$y

  test_res <- test_flexynesis(
    training_result = train_res,
    X_test_list     = X_test,
    y_test          = y_test,
    ...
  )

  # Harrell + Uno, via the shared helper (non-AFT branch).
  metrics <- .compute_test_metrics(
    train_res = list(pred = .flexynesis_pred_vec(train_res$pred, y_train)),
    test_res  = list(pred = .flexynesis_pred_vec(test_res$pred,  y_test)),
    y_train   = y_train,
    y_test    = y_test,
    joint_method = "flexynesis"
  )

  md <- MultiAssayExperiment::metadata(mae_pred)
  if (is.null(md$analysis))            md$analysis <- list()
  if (is.null(md$analysis$flexynesis)) md$analysis$flexynesis <- list()
  md$analysis$flexynesis[[tag]] <- list(
    test    = test_res,
    metrics = metrics
  )
  MultiAssayExperiment::metadata(mae_pred) <- md

  mae_pred
}
