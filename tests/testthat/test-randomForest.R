test_that("fit_model and predict work with randomForest", {
  skip_if_not_installed("randomForest")

  # Use penguins dataset
  library(palmerpenguins)
  penguins <- na.omit(penguins[, -c(2, 7, 8)])
  train_data <- penguins

  # Fit model
  model <- fit_model(train_data, species ~ ., classifier = randomForest::randomForest)

  expect_s3_class(model, "classbound")
  expect_s3_class(model$fit, "randomForest")

  # Predict
  preds <- predict(model, train_data)

  expect_type(preds, "list")
  expect_true("class" %in% names(preds))
  expect_true("probs" %in% names(preds))

  expect_s3_class(preds$class, "factor")
  expect_true(is.matrix(preds$probs))

  expect_length(preds$class, nrow(train_data))
  expect_equal(nrow(preds$probs), nrow(train_data))
})
