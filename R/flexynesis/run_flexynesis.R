run_flexynesis <- function(
    data_path,
    omic_list,
    surv_time,
    surv_event,
    model_class = "DirectPred",
    threads = 10,
    device = "cpu",
    use_cv = TRUE,
    outdir
) {

  conda_bin <- reticulate::conda_binary()

  cmd_args <- c(
    "run", "-n", "flexynesisenv",
    "flexynesis",
    "--data_path", data_path,
    "--model_class", model_class,
    "--surv_time_var", surv_time,
    "--surv_event_var", surv_event,
    "--data_types", paste(omic_list, collapse = ","),
    "--threads", as.character(threads),
    "--device", device,
    "--outdir", outdir
  )

  if (use_cv) cmd_args <- c(cmd_args, "--use_cv")

  message("Running flexynesis command:")
  message(paste(conda_bin, paste(cmd_args, collapse = " ")))

  res <- system2(
    command = conda_bin,
    args = cmd_args,
    stdout = TRUE,
    stderr = TRUE
  )

  # output python
  cat(res, sep = "\n")

  # check output
  csv_path <- file.path(outdir, "job.predicted_labels.csv")
  if (!file.exists(csv_path)) {
    stop(
      "Flexynesis finished but output file not found:\n",
      csv_path,
      "\nCheck Python output above."
    )
  }

  invisible(TRUE)
}
