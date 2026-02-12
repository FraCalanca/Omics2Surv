export_flexynesis_data <- function(
    X_train_list,
    X_test_list,
    y_train,
    y_test,
    outdir
) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(outdir, "train"))
  dir.create(file.path(outdir, "test"))

  write.csv(y_train, file.path(outdir, "train", "clin.csv"))
  write.csv(y_test,  file.path(outdir, "test",  "clin.csv"))

  for (omic in names(X_train_list)) {
    write.csv(t(X_train_list[[omic]]),
              file.path(outdir, "train", paste0(omic, ".csv")))
    write.csv(t(X_test_list[[omic]]),
              file.path(outdir, "test", paste0(omic, ".csv")))
  }
}
