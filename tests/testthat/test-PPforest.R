test_that("fit_model and predict_model work with PPforest", {
  skip_if_not_installed("PPforest")

  # Use penguins dataset
  library(palmerpenguins)
  penguins <- na.omit(penguins[, -c(2, 7, 8)])
  train_data <- penguins[, 2:5]
  train_labels <- penguins$species

  # Fit model
  # For tests, we use small m to be fast
  model <- fit_model(train_data, train_labels, method = "PPforest", m = 10)

  expect_s3_class(model, "classbound_model")
  expect_equal(model$method, "PPforest")

  # Predict
  preds <- predict_model(model, train_data)

  expect_type(preds, "list")
  expect_true("class" %in% names(preds))
  expect_true("probs" %in% names(preds))

  expect_s3_class(preds$class, "factor")
  expect_null(preds$probs)

  expect_length(preds$class, nrow(train_data))
})
