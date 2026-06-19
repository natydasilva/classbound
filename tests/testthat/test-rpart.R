test_that("fit_model and predict_model work end-to-end with rpart", {
  # Simple penguins dataset subset
  library(palmerpenguins)
  penguins <- na.omit(penguins[, -c(2, 7, 8)])
  train_data <- penguins[, 2:5]
  train_labels <- penguins$species

  # 1. Fit Phase
  model <- fit_model(
    data = train_data,
    labels = train_labels,
    method = "rpart"
  )

  # Check standard wrapper class
  expect_s3_class(model, "classbound_model")
  expect_equal(model$method, "rpart")

  # Check that underlying model got built
  expect_s3_class(model$model, "rpart")

  # 2. Predict Phase
  preds <- predict_model(model, newdata = train_data)

  # Output structure guarantees
  expect_type(preds, "list")
  expect_named(preds, c("class", "probs"))
  
  # Predictions output expectations
  expect_equal(length(preds$class), nrow(train_data))
  expect_true(is.factor(preds$class))
  
  # Probabilities output expectations
  expect_equal(nrow(preds$probs), nrow(train_data))
  expect_equal(ncol(preds$probs), length(levels(train_labels)))
})
