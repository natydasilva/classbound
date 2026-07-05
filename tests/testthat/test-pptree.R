test_that("fit_model and predict work with pptree", {
  skip_if_not_installed("PPtreeExt")

  # Use penguins dataset
  library(palmerpenguins)
  penguins <- na.omit(penguins[, -c(2, 7, 8)])
  train_data <- penguins

  # Fit model
  model <- fit_model(train_data, species ~ ., classifier = PPtreeExt::PPtreeExtclass)

  expect_s3_class(model, "classbound")
  expect_s3_class(model$fit, "PPtreeExtclass")

  # Predict
  test_data <- train_data[, -which(names(train_data) == "species")]
  preds <- predict(model, test_data)

  expect_type(preds, "list")
  expect_true("class" %in% names(preds))
  expect_true("probs" %in% names(preds))

  expect_s3_class(preds$class, "factor")
  expect_null(preds$probs)

  expect_length(preds$class, nrow(train_data))
})
