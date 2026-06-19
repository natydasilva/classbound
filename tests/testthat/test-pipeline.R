
test_that("fit_model handles invalid inputs and unsupported methods gracefully", {
  library(palmerpenguins)
  penguins <- na.omit(penguins[, -c(2, 7, 8)])
  train_data <- penguins[, 2:5]
  train_labels <- penguins$species

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
  library(palmerpenguins)
  penguins <- na.omit(penguins[, -c(2, 7, 8)])
  raw_rpart <- rpart::rpart(species ~ ., data = penguins)
  
  expect_error(
    predict_model(raw_rpart, newdata = penguins),
    "must be a 'classbound_model'"
  )
})
