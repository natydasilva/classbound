test_that("fit_model and predict work with PPforest", {
  skip_if_not_installed("PPforest")

  # Use penguins dataset
  library(palmerpenguins)
  penguins <- as.data.frame(na.omit(penguins[, -c(2, 7, 8)]))
  train_data <- penguins

  # Fit model
  # For tests, we use small m to be fast
  model <- fit_model(
    data = train_data,
    formula = species ~ .,
    classifier = PPforest::PPforest,
    interface = "custom",
    fit_args = list(data = train_data, y = "species", m = 10, size.tr = 1, size.p = 1, PPmethod = "LDA")
  )

  expect_s3_class(model, "classbound")
  expect_s3_class(model$fit, "PPforest")

  # Predict (only pass predictors to PPforest predict)
  test_data <- train_data[, -which(names(train_data) == "species")]
  preds <- predict(model, test_data)

  expect_type(preds, "list")
  expect_true("class" %in% names(preds))
  expect_true("probs" %in% names(preds))

  expect_s3_class(preds$class, "factor")
  expect_null(preds$probs)

  expect_length(preds$class, nrow(train_data))
})
