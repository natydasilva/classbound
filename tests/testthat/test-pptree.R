test_that("fit_model and predict_model work with pptree", {
  skip_if_not_installed("PPtreeExt")

  # Use penguins dataset
  library(palmerpenguins)
  penguins <- na.omit(penguins[, -c(2, 7, 8)])
  train_data <- penguins[, 2:5]
  train_labels <- penguins$species

  # Fit model
  model <- fit_model(train_data, train_labels, method = "pptree")

  expect_s3_class(model, "classbound_model")
  expect_equal(model$method, "pptree")

  # Predict
  preds <- predict_model(model, train_data)

  expect_type(preds, "list")
  expect_true("class" %in% names(preds))
  expect_true("probs" %in% names(preds))

  expect_s3_class(preds$class, "factor")
  expect_null(preds$probs)

  expect_length(preds$class, nrow(train_data))
})
