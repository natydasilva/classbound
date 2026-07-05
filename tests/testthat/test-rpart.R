test_that("fit_model and predict work end-to-end with rpart", {
  skip_if_not_installed("rpart")
  # Simple penguins dataset subset
  library(palmerpenguins)
  penguins <- na.omit(penguins[, -c(2, 7, 8)])
  train_data <- penguins

  # 1. Fit Phase
  model <- fit_model(
    data = train_data,
    formula = species ~ bill_length_mm + bill_depth_mm + flipper_length_mm + body_mass_g,
    classifier = rpart::rpart
  )

  # Check standard wrapper class
  expect_s3_class(model, "classbound")

  # Check that underlying model got built
  expect_s3_class(model$fit, "rpart")

  # 2. Predict Phase
  preds <- predict(model, newdata = train_data)

  # Output structure guarantees
  expect_type(preds, "list")
  expect_named(preds, c("class", "probs"))

  # Predictions output expectations
  expect_equal(length(preds$class), nrow(train_data))
  expect_true(is.factor(preds$class))

  # Probabilities output expectations
  expect_equal(nrow(preds$probs), nrow(train_data))
  expect_equal(ncol(preds$probs), length(levels(train_data$species)))
})
