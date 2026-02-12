import_flexynesis_results <- function(csv_path) {
  preds <- readr::read_csv(csv_path)

  train <- preds[preds$split == "train", ]
  test  <- preds[preds$split == "test", ]

  list(
    train = data.frame(
      pred = train$predicted_label,
      row.names = train$sample_id
    ),
    test = data.frame(
      pred = test$predicted_label,
      row.names = test$sample_id
    )
  )
}
