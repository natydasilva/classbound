test_that("fit_model and predict_model work end-to-end with rpart", {
  # Simple Iris dataset subset
  train_data <- iris[, 1:4]
  train_labels <- iris$Species

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

test_that("fit_model handles invalid inputs and unsupported methods gracefully", {
  train_data <- iris[, 1:4]
  train_labels <- iris$Species

  # Missing method
  expect_error(
    fit_model(data = train_data, labels = train_labels),
    "Please specify a classification method"
  )

  # Unsupported method
  expect_error(
    fit_model(data = train_data, labels = train_labels, method = "unsupported_magic_model"),
    "not supported"
  )
})

test_that("predict_model correctly enforces structural constraints", {
  # Unwrapped model passed directly
  raw_rpart <- rpart::rpart(Species ~ ., data = iris)
  
  expect_error(
    predict_model(raw_rpart, newdata = iris),
    "must be a 'classbound_model'"
  )
})
