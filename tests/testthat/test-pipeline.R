
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
